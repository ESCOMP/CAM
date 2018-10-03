module steady_state_tei

  use shr_kind_mod,   only : r8 => shr_kind_r8            ! Real kind to declare variables
  use physics_buffer, only : pbuf_get_index, &            !
                             physics_buffer_desc, &       !
                             pbuf_get_field, &            ! Needed to access physics buffer
                             pbuf_set_field
  use physics_types,  only : physics_state                ! Structures containing physics state variables
  use ppgrid,         only : pcols, pver, pverp, begchunk, endchunk
  use cam_abortutils, only : endrun

  implicit none
  
  private 

  !------------------------
  ! PUBLIC: interfaces 
  !------------------------
  public :: steady_state_tei_init
  public :: steady_state_tei_tend
  
  integer :: indxTe = -1
  integer :: indxTi = -1
  integer :: indxQt = -1
  integer :: indxAR = -1
  integer :: indxO1 = -1
  integer :: indxO2 = -1
  
  real(r8) :: op_mass
  
contains

!==============================================================================

  subroutine steady_state_tei_init(pbuf2d)
  
!-----------------------------------------------------------------------
! Time independent initialization for ionosphere simulation.
!-----------------------------------------------------------------------

    use constituents,   only : cnst_get_ind
    use mo_chem_utls,   only : get_spc_ndx                  ! Routine to get index of adv_mass array for short lived species
    use chem_mods,      only : adv_mass                     ! Array holding mass values for short lived species
    use infnan,   only: nan, assignment(=)

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    real(r8) :: nanval

    call cnst_get_ind( 'O2', indxO2, abort=.true. )
    call cnst_get_ind( 'O',  indxO1, abort=.true. )

    !-------------------------------------------------------------------------------------------------------------------
    !  Get physics buffer indices
    !-------------------------------------------------------------------------------------------------------------------
    indxTe = pbuf_get_index( 'TElec' )
    indxTi = pbuf_get_index( 'TIon' )
    indxQt = pbuf_get_index( 'QTeAur' )
    indxAR = pbuf_get_index( 'AurIPRateSum' )
    
    nanval=nan
    call pbuf_set_field(pbuf2d, indxQt, nanval)
    call pbuf_set_field(pbuf2d, indxAR, nanval)
    
    op_mass = adv_mass(get_spc_ndx('Op'))

  end subroutine steady_state_tei_init


!==============================================================================

  subroutine steady_state_tei_tend(state,istate, dse_tend, pbuf)
    use tei_mod,          only : settei
    use physconst,        only : mbarv                       ! Constituent dependent mbar
    use solar_parms_data, only : f107=>solar_parms_f107      ! 10.7 cm solar flux
    use mo_apex,          only : alatm
    use perf_mod,         only : t_startf, t_stopf           ! timing utils
    use ionos_state_mod,  only : ionos_state
    !-------------------------------------------------------------------------------------
    ! Calculate dry static energy and O+ tendency for extended ionosphere simulation 
    !-------------------------------------------------------------------------------------

    !------------------------------Arguments--------------------------------

    type(physics_state), intent(in)    :: state               ! physics state structure
    type(ionos_state),target,intent(in)    :: istate       ! ionosphere state structure
    real(r8),            intent(inout) :: dse_tend(:,:)       ! dry static energy tendency
    type(physics_buffer_desc),pointer  :: pbuf(:)             ! physics buffer

    !---------------------------Local variables-------------------------------

    real(r8),dimension(pver,state%ncol) :: &
         tn,    &   ! neutral temperature (deg K) 
         o2,    &   ! molecular oxygen (mmr) 
         o1,    &   ! atomic oxygen (mmr)
         n2,    &   ! molecular nitrogen (mmr)
         ne,    &   ! electron density (cm3)
         te,    &   ! electron temperature (from previous time step) (K)
         ti,    &   ! ion temperature (from previous time step) (K)
         op,    &   ! O+ number dens (/cm^3) 
         o2p,   &   ! O2+ number dens (/cm^3) 
         nop,   &   ! NO+ number dens (/cm^3) 
         barm,  &   ! mean molecular weight (g/mole)
         qji_ti     ! joule heating from qjoule_ti (used ui,vi) ev/g/s

    real(r8) :: chi(state%ncol)             ! solar zenith angle (radians)
    real(r8) :: qtot(pver,state%ncol)       ! total ionization rate  (s-1 cm-3)
    real(r8) :: pmid(pver,state%ncol)       ! mid-level press (dyne/cm^2) !  10._r8*pmid(:ncol,k) 
    real(r8) :: pint(pverp,state%ncol)      ! interface press (dyne/cm^2)

    real(r8) :: te_out(pver,state%ncol)     ! output electron temperature (deg K) 
    real(r8) :: ti_out(pver,state%ncol)     ! output ion temperature (deg K)
    real(r8) :: qtotal_out(pver,state%ncol) ! heating rate of neutrals
    
    integer :: lchnk                        ! Chunk number 
    integer :: ncol                         ! Number of columns in chunk
    integer :: i, k

    real(r8), parameter :: evergs = 1.602e-12_r8           ! 1 eV = 1.602e-12 ergs (ergs/eV)
    real(r8), parameter :: sToQConv  = 6.24E15_r8          ! Conversion from J/kg/s to ev/g/s
    real(r8), parameter :: n2_mass = 28._r8                ! N2 molecular weight kg/kmol

    real(r8), pointer :: te_ptr(:,:)                       ! Pointer to electron temperature in pbuf (K)
    real(r8), pointer :: ti_ptr(:,:)                       ! Pointer to ion temperature in pbuf (K)
    real(r8), pointer :: qteaur(:)
    real(r8), pointer :: aurIPRateSum(:,:)                 ! Auroral ion production sum for O2+,O+,N2+
                                                           ! (s-1 cm-3 from module mo_aurora)

    call t_startf ('steady_state_tei_tend')

    !----------------------------------------------------------------
    !  Get number of this chunk
    !----------------------------------------------------------------
    lchnk = state%lchnk
    ncol = state%ncol

    call pbuf_get_field(pbuf, indxTe, te_ptr)
    call pbuf_get_field(pbuf, indxTi, ti_ptr)
    call pbuf_get_field(pbuf, indxQt, qteaur) 
    call pbuf_get_field(pbuf, indxAR, aurIPRateSum)

    do i =1,ncol
       chi(i) = acos(istate%cosZenAngR(i))
       qji_ti(1:pver,i) = dse_tend(i,pver:1:-1)*sToQConv*evergs ! J/kg/s -->  ev/g/s --> ergs/s/g
       ! take the part going into heating the neutrals
       qji_ti(1:pver,i) = qji_ti(1:pver,i) * mbarv(i,pver:1:-1,lchnk) / (mbarv(i,pver:1:-1,lchnk)+op_mass)
       qtot(1:pver,i) = istate%sumIonPRates(i,pver:1:-1) + aurIPRateSum(i,pver:1:-1)
       te(1:pver,i) = te_ptr(i,pver:1:-1)
       ti(1:pver,i) = ti_ptr(i,pver:1:-1)
       tn(1:pver,i) = state%t(i,pver:1:-1)
       pmid(1:pver,i) = state%pmid(i,pver:1:-1) *10._r8 ! dynes/cm2
       pint(1:pverp,i) = state%pint(i,pverp:1:-1)*10._r8
       o2(1:pver,i) = state%q(i,pver:1:-1,indxO2)
       o1(1:pver,i) = state%q(i,pver:1:-1,indxO1)
       n2(1:pver,i) = istate%n2_mmr(i,pver:1:-1)
       barm(1:pver,i) = mbarv(i,pver:1:-1,lchnk)
       ne(1:pver,i) = istate%ndensE(i,pver:1:-1)
       op(1:pver,i) = istate%ndensOp(i,pver:1:-1)
       o2p(1:pver,i) = istate%ndensO2p(i,pver:1:-1)
       nop(1:pver,i) = istate%ndensNOp(i,pver:1:-1)
    end do

    call settei(tn,o2,o1,n2,ne,te,ti,op,o2p,nop,barm,qji_ti, &
                f107, chi, qtot, qteaur(:ncol), alatm(:ncol,lchnk), &
                istate%dipmag(:ncol,1), pmid, pint, 1,pver,ncol, &
                te_out, ti_out, qtotal_out )
    
    do i =1,ncol
       te_ptr(i,pver:1:-1) = te_out(1:pver,i)
       ti_ptr(i,pver:1:-1) = ti_out(1:pver,i)
       dse_tend(i,pver:1:-1) = qtotal_out(1:pver,i) / evergs / sToQConv ! ergs/g/s --> eV/g/s --> J/kg/sec
    end do

    call t_stopf ('steady_state_tei_tend')

  end subroutine steady_state_tei_tend

end module steady_state_tei
