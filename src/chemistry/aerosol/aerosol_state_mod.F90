module aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use physconst, only: pi

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
     procedure :: get_amb_species_numdens
     procedure :: get_cld_species_numdens
     procedure :: coated_frac
     procedure :: mass_mean_radius
     procedure :: watact_mfactor
     procedure(aero_hetfrz_size_wght), deferred :: hetfrz_size_wght
  end type aerosol_state

  ! for state fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  real(r8), parameter :: per_cm3 = 1.e-6_r8 ! factor for m-3 to cm-3 conversions
  real(r8), parameter :: per_m3 = 1.e6_r8   ! factor for cm-3 to m-3 conversions
  real(r8), parameter :: kg2mug = 1.e9_r8   ! factor for kg to micrograms (mug) conversions

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

     !------------------------------------------------------------------------------
     ! returns the volume-weighted fractions of aerosol subset `bin_ndx` that can act
     ! as heterogeneous freezing nuclei
     !------------------------------------------------------------------------------
     function aero_hetfrz_size_wght(self, bin_ndx, ncol, nlev) result(wght)
       import :: aerosol_state, r8
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx             ! bin number
       integer, intent(in) :: ncol                ! number of columns
       integer, intent(in) :: nlev                ! number of vertical levels

       real(r8) :: wght(ncol,nlev)

     end function aero_hetfrz_size_wght

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
  ! returns ambient aerosol number density for a given bin number and species type
  !------------------------------------------------------------------------------
  subroutine get_amb_species_numdens(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, numdens)
    use aerosol_properties_mod, only: aerosol_properties
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: numdens(:,:)         ! species number densities (#/cm^3)

    real(r8), pointer :: num(:,:)
    real(r8) :: type_wght(ncol,nlev)
    real(r8) :: size_wght(ncol,nlev)

    size_wght = self%hetfrz_size_wght(bin_ndx, ncol, nlev)

    call self%icenuc_type_wght_base(bin_ndx, ncol, nlev, species_type, aero_props, rho, type_wght)

    call self%get_ambient_num(bin_ndx, num)

    numdens(:ncol,:) = num(:ncol,:)*rho(:ncol,:)*type_wght(:ncol,:)*size_wght(:ncol,:)*per_cm3

  end subroutine get_amb_species_numdens

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number density for a given bin number and species type
  !------------------------------------------------------------------------------
  subroutine get_cld_species_numdens(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, numdens)
    use aerosol_properties_mod, only: aerosol_properties
    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: numdens(:,:)         ! number densities (#/cm^3)

    real(r8), pointer :: num(:,:)
    real(r8) :: type_wght(ncol,nlev)
    real(r8) :: size_wght(ncol,nlev)

    size_wght = self%hetfrz_size_wght(bin_ndx, ncol, nlev)

    call self%icenuc_type_wght_base(bin_ndx, ncol, nlev, species_type, aero_props, rho, type_wght, cloud_borne=.true.)

    call self%get_cldbrne_num(bin_ndx, num)

    numdens(:ncol,:) = num(:ncol,:)*rho(:ncol,:)*type_wght(:ncol,:)*size_wght(:ncol,:)*per_cm3

  end subroutine get_cld_species_numdens

  !------------------------------------------------------------------------------
  ! returns aerosol type weights for a given aerosol type and bin
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght_base(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, wght, cloud_borne)

    use aerosol_properties_mod, only: aerosol_properties

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wght(:,:)            ! type weights
    logical, optional, intent(in) :: cloud_borne  ! if TRUE cloud-borne aerosols are used
                                                  ! otherwise ambient aerosols are used

    real(r8) :: mass(ncol,nlev)
    real(r8) :: totalmass(ncol,nlev)
    real(r8), pointer :: aer_bin(:,:)

    character(len=aero_name_len) :: spectype, sptype
    integer :: ispc
    logical :: cldbrne

    if (present(cloud_borne)) then
       cldbrne = cloud_borne
    else
       cldbrne = .false.
    end if

    wght(:,:) = 0._r8
    totalmass(:,:) = 0._r8
    mass(:,:)   = 0._r8

    if (species_type=='sulfate_strat') then
       sptype = 'sulfate'
    else
       sptype = species_type
    end if

    do ispc = 1, aero_props%nspecies(bin_ndx)

       if (cldbrne) then
          call self%get_cldbrne_mmr(ispc, bin_ndx, aer_bin)
       else
          call self%get_ambient_mmr(ispc, bin_ndx, aer_bin)
       end if
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

  !------------------------------------------------------------------------------
  ! returns the fraction of particle surface area of aerosol subset `bin_ndx` covered
  ! by at least a monolayer of species `species_type` [0-1]
  !------------------------------------------------------------------------------
  function coated_frac(self, bin_ndx,  species_type, ncol, nlev, aero_props, radius) result(frac)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: radius(:,:) ! m

    real(r8) :: frac(ncol,nlev)                              ! coated fraction

    !------------coated variables--------------------
    real(r8), parameter :: n_so4_monolayers_dust = 1.0_r8 ! number of so4(+nh4) monolayers needed to coat a dust particle
    real(r8), parameter :: dr_so4_monolayers_dust = n_so4_monolayers_dust * 4.76e-10_r8
    real(r8) :: vol_shell(ncol,nlev)
    real(r8) :: vol_core(ncol,nlev)
    real(r8) :: alnsg, fac_volsfc
    real(r8) :: tmp1(ncol,nlev), tmp2(ncol,nlev)
    real(r8),pointer :: sulf_mmr(:,:)
    real(r8),pointer :: soa_mmr(:,:)
    real(r8),pointer :: pom_mmr(:,:)
    real(r8),pointer :: aer_mmr(:,:)

    integer :: sulf_ndx
    integer :: soa_ndx
    integer :: pom_ndx
    integer :: species_ndx

    real(r8) :: specdens_so4
    real(r8) :: specdens_pom
    real(r8) :: specdens_soa
    real(r8) :: specdens

    character(len=aero_name_len) :: spectype
    integer :: ispc

    frac = -huge(1._r8)

    sulf_ndx = -1
    pom_ndx = -1
    soa_ndx = -1
    species_ndx = -1

    do ispc = 1, aero_props%nspecies(bin_ndx)
       call aero_props%species_type(bin_ndx, ispc, spectype)

       select case ( trim(spectype) )
       case('sulfate')
          sulf_ndx = ispc
       case('p-organic')
          pom_ndx = ispc
       case('s-organic')
          soa_ndx = ispc
       end select
       if (spectype==species_type) then
          species_ndx = ispc
       end if
    end do

    vol_shell(:ncol,:) = 0._r8

    if (sulf_ndx>0) then
       call aero_props%get(bin_ndx, sulf_ndx, density=specdens_so4)
       call self%get_ambient_mmr(sulf_ndx, bin_ndx, sulf_mmr)
       vol_shell(:ncol,:) = vol_shell(:ncol,:) + sulf_mmr(:ncol,:)/specdens_so4
    end if
    if (pom_ndx>0) then
       call aero_props%get(bin_ndx, pom_ndx, density=specdens_pom)
       call self%get_ambient_mmr(pom_ndx, bin_ndx, pom_mmr)
       vol_shell(:ncol,:) = vol_shell(:ncol,:) + pom_mmr(:ncol,:)*aero_props%pom_equivso4_factor()/specdens_pom
    end if
    if (soa_ndx>0) then
       call aero_props%get(bin_ndx, soa_ndx, density=specdens_soa)
       call self%get_ambient_mmr(soa_ndx, bin_ndx, soa_mmr)
       vol_shell(:ncol,:) = vol_shell(:ncol,:) + soa_mmr(:ncol,:)*aero_props%soa_equivso4_factor()/specdens_soa
    end if

    call aero_props%get(bin_ndx, species_ndx, density=specdens)
    call self%get_ambient_mmr(species_ndx, bin_ndx, aer_mmr)
    vol_core(:ncol,:) = aer_mmr(:ncol,:)/specdens

    alnsg = aero_props%alogsig(bin_ndx)
    fac_volsfc = exp(2.5_r8*alnsg**2)

    tmp1(:ncol,:) = vol_shell(:ncol,:)*(radius(:ncol,:)*2._r8)*fac_volsfc
    tmp2(:ncol,:) = max(6.0_r8*dr_so4_monolayers_dust*vol_core(:ncol,:), 0.0_r8)

    where(tmp1(:ncol,:)>0._r8 .and. tmp2(:ncol,:)>0._r8)
       frac(:ncol,:) = tmp1(:ncol,:)/tmp2(:ncol,:)
    elsewhere
       frac(:ncol,:) = 0.001_r8
    end where

    where(frac(:ncol,:)>1._r8)
       frac(:ncol,:) = 1._r8
    end where
    where(frac(:ncol,:) < 0.001_r8)
       frac(:ncol,:) = 0.001_r8
    end where

  end function coated_frac

  !------------------------------------------------------------------------------
  ! returns the radius [m] of particles in aerosol subset `bin_ndx` assuming all particles are
  ! the same size and only species `species_ndx` contributes to the particle volume
  !------------------------------------------------------------------------------
  function mass_mean_radius(self, bin_ndx, species_ndx, ncol, nlev, aero_props, rho) result(radius)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: species_ndx            ! species number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)

    real(r8) :: radius(ncol,nlev) ! m

    character(len=aero_name_len) :: species_type
    real(r8) :: aer_numdens(ncol,nlev) ! kg/m3
    real(r8) :: aer_massdens(ncol,nlev) ! kg/m3
    real(r8),pointer :: aer_mmr(:,:) ! kg/kg

    real(r8) :: specdens,minrad
    real(r8) :: wght(ncol,nlev)
    integer :: i,k

    wght = self%hetfrz_size_wght(bin_ndx, ncol, nlev)

    call aero_props%species_type(bin_ndx, species_ndx, spectype=species_type)

    call aero_props%get(bin_ndx, species_ndx, density=specdens) ! kg/m3
    call self%get_ambient_mmr(species_ndx, bin_ndx, aer_mmr) ! kg/kg
    call self%get_amb_species_numdens(bin_ndx, ncol, nlev, species_type, aero_props, rho, aer_numdens) ! #/cm3

    aer_massdens(:ncol,:) = aer_mmr(:ncol,:)*rho(:ncol,:)*wght(:ncol,:) ! kg/m3

    minrad = aero_props%min_mass_mean_rad(bin_ndx, species_ndx)

    do k = 1,nlev
       do i = 1,ncol
          if (aer_massdens(i,k)*1.0e-3_r8 > 1.0e-30_r8 .and. aer_numdens(i,k) > 1.0e-3_r8) then
             radius(i,k) = (3._r8/(4*pi*specdens)*aer_massdens(i,k)/(aer_numdens(i,k)*per_m3))**(1._r8/3._r8) ! m
          else
             radius(i,k) = minrad
          end if
       end do
    end do

  end function mass_mean_radius

  !------------------------------------------------------------------------------
  ! calculates water activity mass factor -- density*(1.-(OC+BC)/(OC+BC+SO4)) [mug m-3]
  ! of species `species_type` in subset `bin_ndx`
  !------------------------------------------------------------------------------
  subroutine watact_mfactor(self, bin_ndx,  species_type, ncol, nlev, aero_props, rho, wact_factor)

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wact_factor(:,:)     ! water activity factor -- density*(1.-(OC+BC)/(OC+BC+SO4)) [mug m-3]

    real(r8), pointer :: aer_mmr(:,:)
    real(r8), pointer :: bin_num(:,:)
    real(r8) :: tot2_mmr(ncol,nlev)
    real(r8) :: tot1_mmr(ncol,nlev)
    real(r8) :: aer_numdens(ncol,nlev)
    integer :: ispc
    character(len=aero_name_len) :: spectype

    real(r8) :: awcam(ncol,nlev)          ! mass density [mug m-3]
    real(r8) :: awfacm(ncol,nlev)         ! mass factor ! (OC+BC)/(OC+BC+SO4)

    tot2_mmr = 0.0_r8
    tot1_mmr = 0.0_r8

    if (aero_props%soluble(bin_ndx)) then

       do ispc = 1, aero_props%nspecies(bin_ndx)

          call aero_props%species_type(bin_ndx, ispc, spectype)

          if (trim(spectype)=='black-c' .or. trim(spectype)=='p-organic' .or. trim(spectype)=='s-organic') then
             call self%get_ambient_mmr(ispc, bin_ndx, aer_mmr)
             tot2_mmr(:ncol,:) = tot2_mmr(:ncol,:) + aer_mmr(:ncol,:)
          end if
          if (trim(spectype)=='sulfate') then
             call self%get_ambient_mmr(ispc, bin_ndx, aer_mmr)
             tot1_mmr(:ncol,:) = tot1_mmr(:ncol,:) + aer_mmr(:ncol,:)
          end if
       end do

    end if

    tot1_mmr(:ncol,:) = tot1_mmr(:ncol,:) + tot2_mmr(:ncol,:)

    call self%get_amb_species_numdens(bin_ndx, ncol, nlev, species_type, aero_props, rho, aer_numdens) ! #/cm3
    call self%get_ambient_num(bin_ndx, bin_num) ! #/kg

    where(bin_num(:ncol,:)>0._r8)
       awcam(:ncol,:) = ((aer_numdens(:ncol,:)*per_m3/bin_num(:ncol,:)) * tot1_mmr(:ncol,:)) * kg2mug  ! [mug m-3]
    elsewhere
       awcam(:ncol,:) = 0._r8
    end where

    where(tot1_mmr(:ncol,:)>0)
       awfacm(:ncol,:) = tot2_mmr(:ncol,:) / tot1_mmr(:ncol,:)
    elsewhere
       awfacm(:ncol,:) = 0._r8
    end where

    wact_factor(:ncol,:) = awcam(:ncol,:)*(1._r8-awfacm(:ncol,:))

  end subroutine watact_mfactor

end module aerosol_state_mod
