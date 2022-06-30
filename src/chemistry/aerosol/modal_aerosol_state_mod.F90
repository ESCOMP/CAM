module modal_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use rad_constituents, only: rad_cnst_get_aer_mmr, rad_cnst_get_mode_num
  use physics_buffer, only: physics_buffer_desc
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties
  use cam_abortutils, only: endrun

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
  function constructor(state,pbuf) result(newobj)
    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(modal_aerosol_state), pointer :: newobj

    character(len=*),parameter :: prefix = 'modal_aerosol_state::constructor: '
    integer :: ierr

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       call endrun(prefix//'error allocating newobj')
    end if

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
  ! returns ambient aerosol mass mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr(self, species_ndx, bin_ndx, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    call rad_cnst_get_aer_mmr(0, bin_ndx, species_ndx, 'a', self%state, self%pbuf, mmr)
  end subroutine get_ambient_mmr

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    call rad_cnst_get_aer_mmr(0, bin_ndx, species_ndx, 'c', self%state, self%pbuf, mmr)
  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer   :: num(:,:)    ! number densities

    call rad_cnst_get_mode_num(0, bin_ndx, 'a', self%state, self%pbuf, num)
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), pointer :: num(:,:)

    call rad_cnst_get_mode_num(0, bin_ndx, 'c', self%state, self%pbuf, num)
  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  ! returns interstitial and cloud-borne aerosol states
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: ibin,ispc, indx

    do ibin = 1, aero_props%nbins()
       indx = aero_props%indexer(ibin, 0)
       call self%get_ambient_num(ibin, raer(indx)%fld)
       call self%get_cldbrne_num(ibin, qqcw(indx)%fld)
       do ispc = 1, aero_props%nspecies(ibin)
          indx = aero_props%indexer(ibin, ispc)
          call self%get_ambient_mmr(ispc,ibin, raer(indx)%fld)
          call self%get_cldbrne_mmr(ispc,ibin, qqcw(indx)%fld)
       end do
    end do

  end subroutine get_states

end module modal_aerosol_state_mod
