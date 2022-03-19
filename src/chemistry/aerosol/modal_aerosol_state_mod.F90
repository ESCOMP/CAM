module modal_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use rad_constituents, only: rad_cnst_get_aer_mmr, rad_cnst_get_mode_num
  use physics_buffer, only: physics_buffer_desc
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: modal_aerosol_state

  type, extends(aerosol_state) :: modal_aerosol_state
     private
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
   contains

     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: get_states

     final :: destructor

  end type modal_aerosol_state

  interface modal_aerosol_state
     procedure :: constructor
  end interface modal_aerosol_state

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf,props) result(newobj)
    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    class(aerosol_properties), intent(in) :: props

    type(modal_aerosol_state), pointer :: newobj

    allocate(newobj)

    newobj%state => state
    newobj%pbuf => pbuf

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_state), intent(inout) :: self

    nullify(self%state)
    nullify(self%pbuf)

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr(self,l,m, x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_aer_mmr(0, m, l, 'a', self%state, self%pbuf, x)
  end subroutine get_ambient_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self,l,m, x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_aer_mmr(0, m, l, 'c', self%state, self%pbuf, x)
  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self,m, x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_mode_num(0, m, 'a', self%state, self%pbuf, x)
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self,m, x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_mode_num(0, m, 'c', self%state, self%pbuf, x)
  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m,mm,l

    do m = 1, aero_props%nbins()
       mm = aero_props%indexer(m, 0)
       call self%get_ambient_num(m, raer(mm)%fld)
       call self%get_cldbrne_num(m, qqcw(mm)%fld)
       do l = 1, aero_props%nspecies(m)
          mm = aero_props%indexer(m, l)
          call self%get_ambient_mmr(l,m, raer(mm)%fld)
          call self%get_cldbrne_mmr(l,m, qqcw(mm)%fld)
       end do
    end do

  end subroutine get_states

end module modal_aerosol_state_mod
