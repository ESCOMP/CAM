!--------------------------------------------------------------------------------
! CARMA diagnostics data object
!--------------------------------------------------------------------------------
module carma_diags_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst
  use ppgrid, only: pcols
  use carma_intr, only: MAXCLDAERDIAG, carma_calculate_cloudborne_diagnostics, carma_output_budget_diagnostics, &
       carma_output_cloudborne_diagnostics
  use carma_flags_mod, only: carma_do_package_diags

  use camsrfexch, only: cam_in_t
  use physics_types, only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc

  implicit none

  private

  public :: carma_diags_t

  !------------------------------------------------------------------------------
  ! CARMA diags object
  !------------------------------------------------------------------------------
  type :: carma_diags_t
     private

     ! CARMA diagnostics
     real(r8), allocatable :: aerclddiag(:,:) ! the cloudborne aerosol diags snapshot
     real(r8), allocatable :: old_cflux(:,:)  ! cam_in%clfux from before the timestep_tend

   contains

     procedure :: update
     procedure :: output

     final :: destructor
  end type carma_diags_t

  interface carma_diags_t
     procedure :: constructor
  end interface carma_diags_t


contains

  !------------------------------------------------------------------------------
  ! object constructor allocates memory
  !------------------------------------------------------------------------------
  function constructor() result(newobj)

    type(carma_diags_t), pointer :: newobj

    integer :: ierr

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    if (.not.carma_do_package_diags) return

    allocate(newobj%aerclddiag(pcols,MAXCLDAERDIAG),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%old_cflux(pcols,pcnst),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

  end function constructor

  !------------------------------------------------------------------------------
  ! update the arrays
  !------------------------------------------------------------------------------
  subroutine update(self, cam_in, state, pbuf)
    class(carma_diags_t), intent(inout) :: self

    type(cam_in_t),      intent(in) :: cam_in
    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    if (.not.carma_do_package_diags) return

    self%old_cflux = cam_in%cflx
    call carma_calculate_cloudborne_diagnostics(state, pbuf, self%aerclddiag)

  end subroutine update

  !------------------------------------------------------------------------------
  ! output the carma bugdets to cam history
  !------------------------------------------------------------------------------
  subroutine output(self, state, ptend, cam_in, label, dt, pbuf)
    class(carma_diags_t), intent(in) :: self

    type(physics_state), intent(in) :: state
    type(physics_ptend), intent(in) :: ptend
    type(cam_in_t),      intent(in) :: cam_in
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: dt
    type(physics_buffer_desc), pointer :: pbuf(:)

    if (.not.carma_do_package_diags) return

    call carma_output_budget_diagnostics(state, ptend, self%old_cflux, cam_in%cflx, dt, label)
    call carma_output_cloudborne_diagnostics(state, pbuf, label, dt, self%aerclddiag)

  end subroutine output

  !------------------------------------------------------------------------------
  ! free up memory
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(carma_diags_t), intent(inout) :: self

    if (allocated(self%aerclddiag)) then
       deallocate(self%aerclddiag)
    end if
    if (allocated(self%old_cflux)) then
       deallocate(self%old_cflux)
    end if

  end subroutine destructor

end module carma_diags_mod
