module aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile, only: iulog
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: aerosol_state
  public :: ptr2d_t

  !> aerosol_state defines the interface to the time-varying aerosol state
  !! variables (e.g., mixing ratios, number concentrations). This includes the
  !! aerosol portion of the overall model state.
  !!
  !! Each aerosol package (e.g., MAM, CARMA, etc) must extend the aerosol_state
  !! class to allow access to the state information (transported and not transported)
  !! of the aerosol package. Any package must implement each of the deferred
  !! procedures of the abstract aerosol_state class, may include additional private
  !! data members and type-bound procedures, and may override functions of the
  !! abstract class.
  !!
  !! Please see the modal_aerosol_state module for an example of how the aerosol_state
  !! class can be extended for a specific aerosol package.
  type, abstract :: aerosol_state
   contains
     procedure(aero_get_transported), deferred :: get_transported
     procedure(aero_set_transported), deferred :: set_transported
     procedure(aero_get_state_mmr), deferred :: get_ambient_mmr
     procedure(aero_get_state_mmr), deferred :: get_cldbrne_mmr
     procedure(aero_get_state_num), deferred :: get_ambient_num
     procedure(aero_get_state_num), deferred :: get_cldbrne_num
     procedure(aero_get_states), deferred :: get_states
     procedure :: loadaer
  end type aerosol_state

  ! for state fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  interface

     !------------------------------------------------------------------------
     ! returns aerosol mass mixing ratio for a given species index and bin index
     !------------------------------------------------------------------------
     subroutine aero_get_state_mmr(self, species_ndx, bin_ndx, mmr)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: species_ndx  ! species index
       integer, intent(in) :: bin_ndx      ! bin index
       real(r8), pointer :: mmr(:,:)       ! mass mixing ratios
     end subroutine aero_get_state_mmr

     !------------------------------------------------------------------------
     ! returns aerosol number mixing ratio for a given species index and bin index
     !------------------------------------------------------------------------
     subroutine aero_get_state_num(self, bin_ndx, num)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx     ! bin index
       real(r8), pointer   :: num(:,:)    ! number densities
     end subroutine aero_get_state_num

     !------------------------------------------------------------------------
     ! returns interstitial and cloud-borne aerosol states
     !------------------------------------------------------------------------
     subroutine aero_get_states( self, aero_props, raer, qqcw )
       import :: aerosol_state, aerosol_properties, ptr2d_t

       class(aerosol_state), intent(in) :: self
       class(aerosol_properties), intent(in) :: aero_props ! properties of the aerosol model
       type(ptr2d_t), intent(out) :: raer(:) ! state of interstitial aerosols
       type(ptr2d_t), intent(out) :: qqcw(:) ! state of cloud-borne aerosols

     end subroutine aero_get_states

     !------------------------------------------------------------------------------
     ! sets transported components
     ! This updates the aerosol model state from the host transported aerosol constituents array.
     ! (mass mixing ratios or number mixing ratios)
     !------------------------------------------------------------------------------
     subroutine aero_set_transported( self, transported_array )
       import :: aerosol_state, r8
       class(aerosol_state), intent(inout) :: self
       real(r8), intent(in) :: transported_array(:,:,:)
     end subroutine aero_set_transported

     !------------------------------------------------------------------------------
     ! returns transported components
     ! This updates the transported aerosol constituent array to match the aerosol model state.
     ! (mass mixing ratios or number mixing ratios)
     !------------------------------------------------------------------------------
     subroutine aero_get_transported( self, transported_array )
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       real(r8), intent(out) :: transported_array(:,:,:)
     end subroutine aero_get_transported

  end interface

contains

  !------------------------------------------------------------------------------
  ! returns aerosol number, volume concentrations, and bulk hygroscopicity
  !------------------------------------------------------------------------------
  subroutine loadaer( self, aero_props, istart, istop, k,  m, cs, phase, &
                       naerosol, vaerosol, hygro)

    use aerosol_properties_mod, only: aerosol_properties

    ! input arguments
    class(aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop       ! stop column index
    integer,  intent(in) :: k           ! level index
    integer,  intent(in) :: m           ! mode or bin index
    real(r8), intent(in) :: cs(:,:)     ! air density (kg/m3)
    integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

    ! output arguments
    real(r8), intent(out) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(out) :: vaerosol(:)  ! volume conc (m3/m3)
    real(r8), intent(out) :: hygro(:)     ! bulk hygroscopicity of mode

    ! internal
    real(r8), pointer :: raer(:,:) ! interstitial aerosol mass, number mixing ratios
    real(r8), pointer :: qqcw(:,:) ! cloud-borne aerosol mass, number mixing ratios
    real(r8) :: specdens, spechygro

    real(r8) :: vol(istart:istop) ! aerosol volume mixing ratio
    integer  :: i, l
    !-------------------------------------------------------------------------------

    do i = istart, istop
       vaerosol(i) = 0._r8
       hygro(i)    = 0._r8
    end do

    do l = 1, aero_props%nspecies(m)

       call self%get_ambient_mmr(l,m, raer)
       call self%get_cldbrne_mmr(l,m, qqcw)
       call aero_props%get(m,l, density=specdens, hygro=spechygro)

       if (phase == 3) then
          do i = istart, istop
             vol(i) = max(raer(i,k) + qqcw(i,k), 0._r8)/specdens
          end do
       else if (phase == 2) then
          do i = istart, istop
             vol(i) = max(qqcw(i,k), 0._r8)/specdens
          end do
       else if (phase == 1) then
          do i = istart, istop
             vol(i) = max(raer(i,k), 0._r8)/specdens
          end do
       else
          write(iulog,*)'phase = ',phase,' in aerosol_state::loadaer not recognized'
          call endrun('phase error in aerosol_state::loadaer')
       end if

       do i = istart, istop
          vaerosol(i) = vaerosol(i) + vol(i)
          hygro(i)    = hygro(i) + vol(i)*spechygro
       end do

    end do

    do i = istart, istop
       if (vaerosol(i) > 1.0e-30_r8) then
          hygro(i)    = hygro(i)/(vaerosol(i))
          vaerosol(i) = vaerosol(i)*cs(i,k)
       else
          hygro(i)    = 0.0_r8
          vaerosol(i) = 0.0_r8
       end if
    end do

    ! aerosol number mixing ratios (#/kg)
    call self%get_ambient_num(m, raer)
    call self%get_cldbrne_num(m, qqcw)
    if (phase == 3) then
       do i = istart, istop
          naerosol(i) = (raer(i,k) + qqcw(i,k))*cs(i,k) ! #/kg -> #/m3
       end do
    else if (phase == 2) then
       do i = istart, istop
          naerosol(i) = qqcw(i,k)*cs(i,k)
       end do
    else
       do i = istart, istop
          naerosol(i) = raer(i,k)*cs(i,k)
       end do
    end if

    ! adjust number
    call aero_props%apply_number_limits( naerosol, vaerosol, istart, istop, m )

  end subroutine loadaer

end module aerosol_state_mod
