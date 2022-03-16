module aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile, only: iulog
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: aerosol_state
  public :: ptr2d_t

  type, abstract :: aerosol_state
   contains
     procedure(aero_get_state_mmr), deferred :: get_ambient_mmr
     procedure(aero_get_state_mmr), deferred :: get_cldbrne_mmr
     procedure(aero_get_state_num), deferred :: get_ambient_num
     procedure(aero_get_state_num), deferred :: get_cldbrne_num
     procedure(aero_get_states), deferred :: get_states
     procedure, private :: loadaer1
     procedure, private :: loadaer2
     generic :: loadaer => loadaer1, loadaer2
  end type aerosol_state

  ! for state fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  interface
     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     function aero_get_state_mmr(self, l,m) result(x)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: l
       integer, intent(in) :: m
       real(r8), pointer :: x(:,:)
     end function aero_get_state_mmr

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     function aero_get_state_num(self, m) result(x)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: m
       real(r8), pointer :: x(:,:)
     end function aero_get_state_num

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_get_states( self, aero_props, raer, qqcw )
       import

       class(aerosol_state), intent(in) :: self
       class(aerosol_properties), intent(in) :: aero_props
       type(ptr2d_t), intent(out) :: raer(:)
       type(ptr2d_t), intent(out) :: qqcw(:)

     end subroutine aero_get_states

  end interface

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine loadaer1( self, aero_props, istart, istop, k,  m, cs, phase, &
                       naerosol, vaerosol, hygro)

    use aerosol_properties_mod, only: aerosol_properties
    use modal_aerosol_properties_mod, only: modal_aerosol_properties

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    ! input arguments
    class(aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop       ! stop column index
    integer,  intent(in) :: m           ! mode or bin index
    integer,  intent(in) :: k           ! level index
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

       raer => self%get_ambient_mmr(l,m)
       qqcw => self%get_cldbrne_mmr(l,m)
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
          write(iulog,*)'phase = ',phase,' in loadaer not recognized'
          call endrun('phase error in loadaer')
       end if

       do i = istart, istop
          vaerosol(i) = vaerosol(i) + vol(i)
          hygro(i)    = hygro(i) + vol(i)*spechygro
       end do

    end do

    do i = istart, istop
       if (vaerosol(i) > 1.0e-30_r8) then   ! +++xl add 8/2/2007
          hygro(i)    = hygro(i)/(vaerosol(i))
          vaerosol(i) = vaerosol(i)*cs(i,k)
       else
          hygro(i)    = 0.0_r8
          vaerosol(i) = 0.0_r8
       end if
    end do

    ! aerosol number
    raer => self%get_ambient_num(m)
    qqcw => self%get_cldbrne_num(m)
    if (phase == 3) then
       do i = istart, istop
          naerosol(i) = (raer(i,k) + qqcw(i,k))*cs(i,k)
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

    select type(aero_props)
    type is (modal_aerosol_properties)

       ! adjust number so that dgnumlo < dgnum < dgnumhi not done for bins
       do i = istart, istop
          naerosol(i) = max(naerosol(i), vaerosol(i)*aero_props%voltonumbhi(m))
          naerosol(i) = min(naerosol(i), vaerosol(i)*aero_props%voltonumblo(m))
       end do
    end select

  end subroutine loadaer1

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine loadaer2( self, aero_props, i, k, m, cs, phase, &
                       naerosol, vaerosol, hygro )

    use aerosol_properties_mod, only: aerosol_properties
    use ppgrid, only: pcols, pver

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    ! input arguments
    class(aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer,  intent(in) :: i           ! column index
    integer,  intent(in) :: k           ! level index
    integer,  intent(in) :: m           ! mode or bin index
    real(r8), intent(in) :: cs          ! air density (kg/m3)
    integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

    ! output arguments
    real(r8), intent(out) :: naerosol  ! number conc (1/m3)
    real(r8), intent(out) :: vaerosol  ! volume conc (m3/m3)
    real(r8), intent(out) :: hygro     ! bulk hygroscopicity of mode

    real(r8) :: cs_a(pcols,pver)          ! air density (kg/m3)
    real(r8) :: naerosol_a(pcols)  ! number conc (1/m3)
    real(r8) :: vaerosol_a(pcols)  ! volume conc (m3/m3)
    real(r8) :: hygro_a(pcols)     ! bulk hygroscopicity of mode

    cs_a(i,k) = cs

    call self%loadaer(aero_props, i, i, k, m, cs_a, phase, naerosol_a, vaerosol_a, hygro_a)

    naerosol = naerosol_a(i)
    vaerosol = vaerosol_a(i)
    hygro = hygro_a(i)

  end subroutine loadaer2

end module aerosol_state_mod
