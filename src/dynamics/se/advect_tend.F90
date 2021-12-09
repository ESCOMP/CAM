!----------------------------------------------------------------------
! this module computes the total advection tendencies of advected
! constituents for the finite volume dycore
!----------------------------------------------------------------------
module advect_tend

  use shr_kind_mod, only : r8 => shr_kind_r8

  save
  private

  public :: compute_adv_tends_xyz

  real(r8), allocatable :: adv_tendxyz(:,:,:,:,:)

contains

  !----------------------------------------------------------------------
  ! computes the total advective tendencies
  ! called twice each time step:
  !   - first call sets the initial mixing ratios
  !   - second call computes and outputs the tendencies
  !----------------------------------------------------------------------
  subroutine compute_adv_tends_xyz(elem,fvm,nets,nete,qn0,n0)
    use cam_history,            only: outfld, hist_fld_active
    use time_manager,           only: get_step_size
    use constituents,           only: tottnam,pcnst    
    use dimensions_mod,         only: nc,np,nlev,ntrac
    use element_mod,            only: element_t
    use fvm_control_volume_mod, only: fvm_struct    
    implicit none

    type (element_t), intent(in) :: elem(:)
    type(fvm_struct), intent(in) :: fvm(:)
    integer,          intent(in) :: nets,nete,qn0,n0
    real(r8) :: dt,idt
    integer  :: i,j,ic,nx,ie
    logical  :: init
    real(r8), allocatable, dimension(:,:) :: ftmp

    if (ntrac>0) then
      nx=nc
    else
      nx=np
    endif
    allocate( ftmp(nx*nx,nlev) )
    
    init = .false.
    if ( .not. allocated( adv_tendxyz ) ) then
      init = .true.
      allocate( adv_tendxyz(nx,nx,nlev,pcnst,nets:nete) )
      adv_tendxyz(:,:,:,:,:) = 0._r8
    endif

    if (ntrac>0) then
      do ie=nets,nete
        do ic=1,pcnst
          adv_tendxyz(:,:,:,ic,ie) = fvm(ie)%c(1:nc,1:nc,:,ic) - adv_tendxyz(:,:,:,ic,ie)
        end do
      end do
    else
      do ie=nets,nete
        do ic=1,pcnst
          adv_tendxyz(:,:,:,ic,ie) = elem(ie)%state%Qdp(:,:,:,ic,qn0)/elem(ie)%state%dp3d(:,:,:,n0)  - adv_tendxyz(:,:,:,ic,ie)
        enddo
      end do
    end if

    if ( .not. init ) then
      dt = get_step_size()
      idt = 1._r8/dt

      do ie=nets,nete
        do ic = 1,pcnst
          do j=1,nx
            do i=1,nx
              ftmp(i+(j-1)*nx,:) = adv_tendxyz(i,j,:,ic,ie)
            end do
          end do
          call outfld(tottnam(ic), ftmp,nx*nx, ie)
        end do
      end do
      deallocate(adv_tendxyz)
    endif
    deallocate(ftmp)
  end subroutine compute_adv_tends_xyz

end module advect_tend
