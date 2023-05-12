module aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len

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
     procedure(aero_get_amb_total_bin_mmr), deferred :: ambient_total_bin_mmr
     procedure(aero_get_state_mmr), deferred :: get_ambient_mmr
     procedure(aero_get_state_mmr), deferred :: get_cldbrne_mmr
     procedure(aero_get_state_num), deferred :: get_ambient_num
     procedure(aero_get_state_num), deferred :: get_cldbrne_num
     procedure(aero_get_states), deferred :: get_states
     procedure(aero_update_bin), deferred :: update_bin
     procedure :: loadaer
     procedure(aero_icenuc_size_wght_arr), deferred :: icenuc_size_wght_arr
     procedure(aero_icenuc_size_wght_val), deferred :: icenuc_size_wght_val
     generic :: icenuc_size_wght => icenuc_size_wght_arr,icenuc_size_wght_val
     procedure :: icenuc_type_wght_base
     procedure :: icenuc_type_wght => icenuc_type_wght_base
     procedure :: nuclice_get_numdens
  end type aerosol_state

  ! for state fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  abstract interface

     !------------------------------------------------------------------------
     ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
     !------------------------------------------------------------------------
     function aero_get_amb_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
       import :: aerosol_state, aerosol_properties, r8
       class(aerosol_state), intent(in) :: self
       class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
       integer, intent(in) :: bin_ndx      ! bin index
       integer, intent(in) :: col_ndx      ! column index
       integer, intent(in) :: lyr_ndx      ! vertical layer index

       real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species

     end function aero_get_amb_total_bin_mmr

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

     !------------------------------------------------------------------------------
     ! return aerosol bin size weights for a given bin
     !------------------------------------------------------------------------------
     subroutine aero_icenuc_size_wght_arr(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx             ! bin number
       integer, intent(in) :: ncol                ! number of columns
       integer, intent(in) :: nlev                ! number of vertical levels
       character(len=*), intent(in) :: species_type  ! species type
       logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
       real(r8), intent(out) :: wght(:,:)

     end subroutine aero_icenuc_size_wght_arr

     !------------------------------------------------------------------------------
     ! return aerosol bin size weights for a given bin, column and vertical layer
     !------------------------------------------------------------------------------
     subroutine aero_icenuc_size_wght_val(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx                ! bin number
       integer, intent(in) :: col_ndx                ! column index
       integer, intent(in) :: lyr_ndx                ! vertical layer index
       character(len=*), intent(in) :: species_type  ! species type
       logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
       real(r8), intent(out) :: wght

     end subroutine aero_icenuc_size_wght_val

     !------------------------------------------------------------------------------
     ! updates state and tendency
     !------------------------------------------------------------------------------
     subroutine aero_update_bin( self, bin_ndx, col_ndx, lyr_ndx, delmmr_sum, delnum_sum, tnd_ndx, dtime, tend )
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx                ! bin number
       integer, intent(in) :: col_ndx                ! column index
       integer, intent(in) :: lyr_ndx                ! vertical layer index
       real(r8),intent(in) :: delmmr_sum             ! mass mixing ratio change summed over all species in bin
       real(r8),intent(in) :: delnum_sum             ! number mixing ratio change summed over all species in bin
       integer, intent(in) :: tnd_ndx                ! tendency index
       real(r8),intent(in) :: dtime                  ! time step size (sec)
       real(r8),intent(inout) :: tend(:,:,:)         ! tendency

     end subroutine aero_update_bin

  end interface

contains

  !------------------------------------------------------------------------------
  ! returns aerosol number, volume concentrations, and bulk hygroscopicity
  !------------------------------------------------------------------------------
  subroutine loadaer( self, aero_props, istart, istop, k,  m, cs, phase, &
                       naerosol, vaerosol, hygro, errnum, errstr)

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

    integer ,         intent(out) :: errnum
    character(len=*), intent(out) :: errstr

    ! internal
    real(r8), pointer :: raer(:,:) ! interstitial aerosol mass, number mixing ratios
    real(r8), pointer :: qqcw(:,:) ! cloud-borne aerosol mass, number mixing ratios
    real(r8) :: specdens, spechygro

    real(r8) :: vol(istart:istop) ! aerosol volume mixing ratio
    integer  :: i, l
    !-------------------------------------------------------------------------------
    errnum = 0

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
          errnum = -1
          write(errstr,*)'phase = ',phase,' in aerosol_state::loadaer not recognized'
          return
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

  !------------------------------------------------------------------------------
  ! returns aerosol type weights for a given aerosol type and bin
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght_base(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, wght)

    use aerosol_properties_mod, only: aerosol_properties

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)        ! air density (kg m-3)
    real(r8), intent(out) :: wght(:,:)            ! type weights

    real(r8) :: mass(ncol,nlev)
    real(r8) :: totalmass(ncol,nlev)
    real(r8), pointer :: aer_bin(:,:)

    character(len=aero_name_len) :: spectype, sptype
    integer :: ispc

    wght(:,:) = 0._r8
    totalmass(:,:) = 0._r8
    mass(:,:)   = 0._r8

    if (species_type=='sulfate_strat') then
       sptype = 'sulfate'
    else
       sptype = species_type
    end if

    do ispc = 1, aero_props%nspecies(bin_ndx)

       call self%get_ambient_mmr(ispc, bin_ndx, aer_bin)
       call aero_props%species_type(bin_ndx, ispc, spectype=spectype)

       totalmass(:ncol,:) = totalmass(:ncol,:) + aer_bin(:ncol,:)*rho(:ncol,:)

       if (trim(spectype) == trim(sptype)) then
          mass(:ncol,:) = mass(:ncol,:) + aer_bin(:ncol,:)*rho(:ncol,:)
       end if

    end do

    where (totalmass(:ncol,:) > 0._r8)
       wght(:ncol,:) = mass(:ncol,:)/totalmass(:ncol,:)
    end where

  end subroutine icenuc_type_wght_base

  !------------------------------------------------------------------------------
  subroutine nuclice_get_numdens(self, aero_props, use_preexisting_ice, ncol, nlev, rho, dust_num_col, sulf_num_col, soot_num_col, sulf_num_tot_col )

    class(aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

    logical, intent(in) :: use_preexisting_ice
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    real(r8), intent(in) :: rho(:,:) ! air density (kg m-3)
    real(r8), intent(out) :: dust_num_col(:,:) ! dust number densities (#/cm^3)
    real(r8), intent(out) :: sulf_num_col(:,:) ! sulfate number densities (#/cm^3)
    real(r8), intent(out) :: soot_num_col(:,:) ! soot number densities (#/cm^3)
    real(r8), intent(out) :: sulf_num_tot_col(:,:) ! stratopsheric sulfate number densities (#/cm^3)

    integer :: ibin,ispc
    character(len=aero_name_len) :: spectype
    real(r8) :: size_wghts(ncol,nlev)
    real(r8) :: type_wghts(ncol,nlev)

    real(r8), pointer :: num_col(:,:)

    real(r8), parameter :: per_cm3 = 1.e-6_r8 ! factor for m-3 to cm-3 conversions

    dust_num_col(:,:) = 0._r8
    sulf_num_col(:,:) = 0._r8
    soot_num_col(:,:) = 0._r8
    sulf_num_tot_col(:,:) = 0._r8

    ! collect number densities (#/cm^3) for dust, sulfate, and soot
    do ibin = 1,aero_props%nbins()

       call self%get_ambient_num(ibin, num_col)

       do ispc = 1,aero_props%nspecies(ibin)

          call aero_props%species_type(ibin, ispc, spectype)

          call self%icenuc_size_wght(ibin, ncol, nlev, spectype, use_preexisting_ice, size_wghts)

          call self%icenuc_type_wght(ibin, ncol, nlev, spectype, aero_props, rho, type_wghts)

          select case ( trim(spectype) )
          case('dust')
             dust_num_col(:ncol,:) = dust_num_col(:ncol,:) &
                  + size_wghts(:ncol,:)*type_wghts(:ncol,:)*num_col(:ncol,:)*rho(:ncol,:)*per_cm3
          case('sulfate')
             ! This order of ops gives bit-for-bit results for cam5 phys ( use_preexisting_ice = .false. )
             sulf_num_col(:ncol,:) = sulf_num_col(:ncol,:) &
                  + num_col(:ncol,:)*rho(:ncol,:)*per_cm3  * size_wghts(:ncol,:)*type_wghts(:ncol,:)
          case('black-c')
             soot_num_col(:ncol,:) = soot_num_col(:ncol,:) &
                  + size_wghts(:ncol,:)*type_wghts(:ncol,:)*num_col(:ncol,:)*rho(:ncol,:)*per_cm3
          end select

       enddo

       ! stratospheric sulfates -- special case not included in the species loop above
       call self%icenuc_size_wght(ibin, ncol, nlev, 'sulfate_strat', use_preexisting_ice, size_wghts)
       call self%icenuc_type_wght(ibin, ncol, nlev, 'sulfate_strat', aero_props, rho, type_wghts)
       sulf_num_tot_col(:ncol,:) = sulf_num_tot_col(:ncol,:) &
            + size_wghts(:ncol,:)*type_wghts(:ncol,:)*num_col(:ncol,:)*rho(:ncol,:)*per_cm3

    enddo

  end subroutine nuclice_get_numdens

end module aerosol_state_mod
