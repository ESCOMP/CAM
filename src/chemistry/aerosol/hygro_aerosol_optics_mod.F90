!-------------------------------------------------------------------------------
! Short-wave hygroscopic aerosol, Long-wave non-hygroscopic (insoluble)
! aerosol optical properties
!-------------------------------------------------------------------------------
module hygro_aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  use aerosol_optics_mod, only: aerosol_optics
  use aerosol_properties_mod, only: aerosol_properties
  use aerosol_state_mod, only: aerosol_state

  implicit none

  private

  public :: hygro_aerosol_optics

  type, extends(aerosol_optics) :: hygro_aerosol_optics

     ! aerosol optics properties tables (from physprops files)
     real(r8), pointer :: ext_sw(:,:) => null()
     real(r8), pointer :: ssa_sw(:,:) => null()
     real(r8), pointer :: asm_sw(:,:) => null()
     real(r8), pointer :: abs_lw(:) => null()

     ! from state
     real(r8), allocatable :: wrh(:,:) ! (-) weighting on left side values ! (pcols,pver)
     integer , allocatable :: krh(:,:) ! index into rh mesh

     ! aerosol mass mixing ratio
     real(r8), pointer :: mmr(:,:)

   contains

     procedure :: sw_props
     procedure :: lw_props

     final :: destructor

  end type hygro_aerosol_optics

  interface hygro_aerosol_optics
     procedure :: constructor
  end interface hygro_aerosol_optics

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(aero_props, aero_state, ibin, ncols, nlevs, numrh, relhum) &
                result(newobj)
    class(aerosol_properties),intent(in) :: aero_props   ! aerosol_properties object
    class(aerosol_state),     intent(in) :: aero_state      ! aerosol_state object
    integer, intent(in) :: ibin   ! bin number
    integer, intent(in) :: ncols, nlevs, numrh
    real(r8),intent(in) :: relhum(ncols,nlevs)

    type(hygro_aerosol_optics), pointer :: newobj

    real(r8) :: rhtrunc(ncols,nlevs)
    integer :: ierr

    allocate(newobj, stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%wrh(ncols,nlevs), stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%krh(ncols,nlevs), stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

! NOTE should try to use table_interp_mod utility !!!
    rhtrunc(1:ncols,1:nlevs)  = min(relhum(1:ncols,1:nlevs),1._r8)
    newobj%krh(1:ncols,1:nlevs) = min(floor( rhtrunc(1:ncols,1:nlevs) * numrh ) + 1, numrh - 1)  ! index into rh mesh
    newobj%wrh(1:ncols,1:nlevs) = rhtrunc(1:ncols,1:nlevs) * numrh - newobj%krh(1:ncols,1:nlevs) ! (-) weighting on left side values

    ! optical properties tables
    call aero_props%optics_params(ibin, &
         sw_hygroscopic_ext=newobj%ext_sw, &
         sw_hygroscopic_ssa=newobj%ssa_sw, &
         sw_hygroscopic_asm=newobj%asm_sw, &
         lw_insoluble_ext=newobj%abs_lw )

    call aero_state%get_ambient_mmr(species_ndx=1, bin_ndx=ibin, mmr=newobj%mmr)

  end function constructor

  !------------------------------------------------------------------------------
  ! returns short wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)

    class(hygro_aerosol_optics), intent(in) :: self

    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pext(ncol) ! parameterized specific extinction (m2/kg)
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)
    real(r8),intent(out) :: palb(ncol) ! parameterized single scattering albedo
    real(r8),intent(out) :: pasm(ncol) ! parameterized asymmetry factor

    integer :: icol

    ! interpolate the properties tables
    do icol = 1, ncol
       pext(icol) = (1._r8 + self%wrh(icol,ilev)) * self%ext_sw(self%krh(icol,ilev)+1,iwav) &
                           - self%wrh(icol,ilev)  * self%ext_sw(self%krh(icol,ilev),  iwav)
       palb(icol) = (1._r8 + self%wrh(icol,ilev)) * self%ssa_sw(self%krh(icol,ilev)+1,iwav) &
                           - self%wrh(icol,ilev)  * self%ssa_sw(self%krh(icol,ilev),  iwav)
       pasm(icol) = (1._r8 + self%wrh(icol,ilev)) * self%asm_sw(self%krh(icol,ilev)+1,iwav) &
                           - self%wrh(icol,ilev)  * self%asm_sw(self%krh(icol,ilev),  iwav)

       pext(icol) = pext(icol) * self%mmr(icol,ilev)

       pabs(icol) = pext(icol) * ( 1._r8 - palb(icol) )

    end do

  end subroutine sw_props

  !------------------------------------------------------------------------------
  ! returns long wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine lw_props(self, ncol, ilev, iwav, pabs)

    class(hygro_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)

    integer :: icol

    pabs(:ncol) = self%abs_lw(iwav) * self%mmr(:ncol,ilev)

  end subroutine lw_props

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)

    type(hygro_aerosol_optics), intent(inout) :: self

    deallocate(self%wrh)
    deallocate(self%krh)

  end subroutine destructor

end module hygro_aerosol_optics_mod
