module bulk_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  !REMOVECAM
  use aerosol_mmr_cam, only: rad_cnst_get_aer_mmr
  !REMOVECAM_END
  !REMOVECAM: no longer need pbuf and state after CAM is retired
  use physics_buffer, only: physics_buffer_desc
  use physics_types,  only: physics_state
  !REMOVECAM_END
  use cam_abortutils,   only: endrun
  use radiative_aerosol, only: rad_aer_get_props
  use string_utils, only: to_lower

  use aerosol_state_mod,      only: aerosol_state, ptr2d_t
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  ! BAM sulfate scaling factor:
  ! former microp_aero_bulk_scale namelist parameter (always 2.0).
  real(r8), parameter :: bam_sulfate_scale = 2.0_r8

  public :: bulk_aerosol_state

  type, extends(aerosol_state) :: bulk_aerosol_state
     private

      !REMOVECAM: state and pbuf will be replaced by SIMA MMR API
      type(physics_state), pointer :: state => null()
      type(physics_buffer_desc), pointer :: pbuf(:) => null()
      !REMOVECAM_END

      ! Per-object workspace for derived number mixing ratio.
      ! Allocated in constructor, deallocated in destructor.
      real(r8), pointer :: num_work_(:,:) => null()   ! (horizontal_dimension, vertical_layer_dimension)
      real(r8), pointer :: zero_fld_(:,:) => null()   ! (horizontal_dimension, vertical_layer_dimension)

   contains

     procedure :: get_transported
     procedure :: set_transported
     procedure :: ambient_total_bin_mmr
     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: get_states
     procedure :: icenuc_size_wght_arr
     procedure :: icenuc_size_wght_val
     procedure :: icenuc_type_wght
     procedure :: update_bin
     procedure :: hetfrz_size_wght
     procedure :: hygroscopicity
     procedure :: water_uptake
     procedure :: dry_volume
     procedure :: wet_volume
     procedure :: water_volume
     procedure :: wet_diameter
     procedure :: convcld_actfrac
     procedure :: wgtpct
     procedure :: get_bulk_num_and_mass
     ! for bit-for-bit
     procedure :: nuclice_get_numdens => nuclice_get_numdens_bam

     final :: destructor

  end type bulk_aerosol_state

  interface bulk_aerosol_state
     procedure :: constructor
  end interface bulk_aerosol_state

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(ncol,state,pbuf,list_idx) result(newobj)
    !REMOVECAM: host-model specific dimensions
    use ppgrid,           only: pcols, pver
    !REMOVECAM_END

    integer, intent(in) :: ncol
    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer, intent(in), optional :: list_idx
    type(bulk_aerosol_state), pointer :: newobj

    integer :: ierr

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    newobj%state => state
    newobj%pbuf => pbuf

    ! set number of active columns internally to prevent loops from accessing beyond
    ! meaningful data in arrays
    call newobj%set_ncol(ncol)

    if (present(list_idx)) call newobj%set_list_idx(list_idx)

    ! Allocate per-object workspace for derived number fields.
    ! Thread-safe: in CAM, each chunk has its own state object.
    allocate(newobj%num_work_(pcols, pver), stat=ierr)
    if (ierr /= 0) call endrun('bulk_aerosol_state constructor: num_work_ allocation error')
    newobj%num_work_(:,:) = 0._r8
    allocate(newobj%zero_fld_(pcols, pver), stat=ierr)
    if (ierr /= 0) call endrun('bulk_aerosol_state constructor: zero_fld_ allocation error')
    newobj%zero_fld_(:,:) = 0._r8

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(bulk_aerosol_state), intent(inout) :: self

    nullify(self%state)
    nullify(self%pbuf)

    if (associated(self%num_work_)) then
       deallocate(self%num_work_)
       nullify(self%num_work_)
    end if
    if (associated(self%zero_fld_)) then
       deallocate(self%zero_fld_)
       nullify(self%zero_fld_)
    end if

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! sets transported components
  ! This aerosol model with the state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(bulk_aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  ! This returns to current state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(bulk_aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine get_transported

  !------------------------------------------------------------------------
  ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
  !------------------------------------------------------------------------
  function ambient_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: bin_ndx      ! bin index
    integer, intent(in) :: col_ndx      ! column index
    integer, intent(in) :: lyr_ndx      ! vertical layer index

    real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    call self%get_ambient_mmr(species_ndx=1, bin_ndx=bin_ndx, mmr=mmr)

    mmr_tot = mmr(col_ndx, lyr_ndx)

  end function ambient_total_bin_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr(self, species_ndx, bin_ndx, mmr)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    ! species_ndx is ignored in the bulk implementation.
    ! bin_ndx is used to identify each individual bulk aerosol.
    call rad_cnst_get_aer_mmr(self%list_idx_, bin_ndx, self%state, self%pbuf, mmr)

  end subroutine get_ambient_mmr

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    ! BAM has no cloud-borne aerosol equivalent, return zero array.
    mmr => self%zero_fld_

  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer   :: num(:,:)    ! number mixing ratio (#/kg)

    real(r8), pointer :: mmr(:,:)
    real(r8)          :: ntm
    character(len=32) :: aname
    integer           :: nc

    ! Derive number mixing ratio from mass: num = mmr * num_to_mass_aer (* bam_sulfate_scale for sulfate).
    ! This matches the inline computation formerly in microp_aero.F90 and nucleate_ice_cam.F90.
    ! Computed into per-object workspace (num_work_); callers must use or copy before the next call.
    ! Only active columns (1:ncol) are computed to avoid FPE on uninitialised padding columns.

    nc = self%ncol()

    call self%get_ambient_mmr(species_ndx=1, bin_ndx=bin_ndx, mmr=mmr)
    call rad_aer_get_props(self%list_idx_, bin_ndx, num_to_mass_aer=ntm, aername=aname)

    ! Apply bam_sulfate_scale to sulfate/volcanic aerosol
    select case ( to_lower( aname(:4) ) )
    case ('sulf', 'volc') ! both treated as 'sulfate' in aero_props%get type.
       self%num_work_(:nc,:) = mmr(:nc,:) * ntm * bam_sulfate_scale
    case default
       self%num_work_(:nc,:) = mmr(:nc,:) * ntm
    end select

    num => self%num_work_

  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), pointer :: num(:,:)

    ! BAM has no cloud-borne equivalent, return zero array.
    num => self%zero_fld_

  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  ! returns interstitial and cloud-borne aerosol states
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    call endrun('ERROR: bulk_aerosol_state_mod%get_states not yet implemented')

  end subroutine get_states

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_arr(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
    real(r8), intent(out) :: wght(:,:)

    ! Empirical 1/25 scaling factor for BAM ice nucleation number densities.
    ! This was previously hardcoded inline in nucleate_ice_cam.F90:633.
    wght(:ncol,:nlev) = 1._r8 / 25._r8

  end subroutine icenuc_size_wght_arr

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin, column and vertical layer
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_val(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
    real(r8), intent(out) :: wght

    ! Empirical 1/25 scaling factor for BAM ice nucleation number densities.
    wght = 1._r8 / 25._r8

  end subroutine icenuc_size_wght_val

  !------------------------------------------------------------------------------
  ! returns aerosol type weights for a given aerosol type and bin
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, wght, cloud_borne)

    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wght(:,:)            ! type weights
    logical, optional, intent(in) :: cloud_borne  ! if TRUE cloud-borne aerosols are used
                                                  ! otherwise ambient aerosols are used

    character(len=32) :: bin_spectype

    ! BAM has exactly 1 species per bin. The type weight is 1.0 when the queried
    ! species type matches the bin's species, 0.0 otherwise. This avoids the
    ! base class computation (which reads MMR just to compute mass/totalmass = 1.0).

    call aero_props%species_type(bin_ndx, 1, bin_spectype)

    if (trim(bin_spectype) == trim(species_type) .or. &
        (species_type == 'sulfate_strat' .and. bin_spectype == 'sulfate')) then
       wght(:ncol,:nlev) = 1._r8
    else
       wght(:ncol,:nlev) = 0._r8
    end if

  end subroutine icenuc_type_wght

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine update_bin( self, bin_ndx, col_ndx, lyr_ndx, delmmr_sum, delnum_sum, tnd_ndx, dtime, tend )
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    real(r8),intent(in) :: delmmr_sum             ! mass mixing ratio change summed over all species in bin
    real(r8),intent(in) :: delnum_sum             ! number mixing ratio change summed over all species in bin
    integer, intent(in) :: tnd_ndx                ! tendency index
    real(r8),intent(in) :: dtime                  ! time step size (sec)
    real(r8),intent(inout) :: tend(:,:,:)         ! tendency

    ! No-op for BAM: ice nucleation does not produce aerosol tendencies
    ! (no interstitial-to-cloud-borne transfer for bulk aerosols).

  end subroutine update_bin

  !------------------------------------------------------------------------------
  ! returns the volume-weighted fractions of aerosol subset `bin_ndx` that can act
  ! as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_size_wght(self, bin_ndx, ncol, nlev) result(wght)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    real(r8) :: wght(ncol,nlev)

    call endrun('ERROR: bulk_aerosol_state_mod%hetfrz_size_wght not yet implemented')

  end function hetfrz_size_wght

  !------------------------------------------------------------------------------
  ! returns hygroscopicity for a given radiation diagnostic list number and
  ! bin number
  !------------------------------------------------------------------------------
  subroutine hygroscopicity(self, bin_ndx, kappa)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx         ! bin number
    real(r8), intent(out) :: kappa(:,:)    ! hygroscopicity (ncol,nlev)

    call endrun('ERROR: bulk_aerosol_state_mod%hygroscopicity not yet implemented')

  end subroutine hygroscopicity

  !------------------------------------------------------------------------------
  ! returns aerosol wet diameter and aerosol water concentration for a given
  ! radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  subroutine water_uptake(self, aero_props, bin_idx, ncol, nlev, dgnumwet, qaerwat)

    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: bin_idx              ! bin number
    integer, intent(in) :: ncol                 ! number of columns
    integer, intent(in) :: nlev                 ! number of levels
    real(r8),intent(out) :: dgnumwet(ncol,nlev) ! aerosol wet diameter (m)
    real(r8),intent(out) :: qaerwat(ncol,nlev)  ! aerosol water concentration (g/g)

    call endrun('ERROR: bulk_aerosol_state_mod%water_uptake not yet implemented')

  end subroutine water_uptake

  !------------------------------------------------------------------------------
  ! aerosol dry volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function dry_volume(self, aero_props, bin_idx, ncol, nlev) result(vol)

    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg
    real(r8), pointer :: mmr(:,:)    ! kg/kg
    real(r8) :: dens                 ! kg/m3

    call aero_props%get(bin_idx, 1, density=dens)
    call self%get_ambient_mmr(species_ndx=1, bin_ndx=bin_idx, mmr=mmr)

    vol(:ncol,:nlev) = mmr(:ncol,:nlev)/dens

  end function dry_volume

  !------------------------------------------------------------------------------
  ! aerosol wet volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function wet_volume(self, aero_props, bin_idx, ncol, nlev) result(vol)

    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    vol = self%dry_volume(aero_props, bin_idx, ncol, nlev) &
        + self%water_volume(aero_props, bin_idx, ncol, nlev)

  end function wet_volume

  !------------------------------------------------------------------------------
  ! aerosol water volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function water_volume(self, aero_props, bin_idx, ncol, nlev) result(vol)

    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    vol = 0._r8

  end function water_volume

  !------------------------------------------------------------------------------
  ! aerosol wet diameter
  !------------------------------------------------------------------------------
  function wet_diameter(self, bin_idx, ncol, nlev) result(diam)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: diam(ncol,nlev)

    call endrun('ERROR: bulk_aerosol_state_mod%wet_diameter not yet implemented')

  end function wet_diameter

  !------------------------------------------------------------------------------
  ! prescribed aerosol activation fraction for convective cloud
  !------------------------------------------------------------------------------
  function convcld_actfrac(self, aero_props, ibin, ispc, ncol, nlev) result(frac)

    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    integer, intent(in) :: ibin   ! bin index
    integer, intent(in) :: ispc   ! species index
    integer, intent(in) :: ncol   ! number of columns
    integer, intent(in) :: nlev   ! number of vertical levels

    real(r8) :: frac(ncol,nlev)

    call endrun('ERROR: bulk_aerosol_state_mod%convcld_actfrac not yet implemented')

  end function convcld_actfrac

  !------------------------------------------------------------------------------
  ! aerosol weight percent of H2SO4/H2O solution
  !------------------------------------------------------------------------------
  function wgtpct(self, ncol, nlev) result(wtp)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in) ::  ncol, nlev
    real(r8) :: wtp(ncol,nlev)  ! weight percent of H2SO4/H2O solution for given icol, ilev

    wtp = -huge(1._r8)

  end function wgtpct

  !------------------------------------------------------------------------------
  ! Compute BAM number concentration (#/m3) and mass concentration (kg/m3)
  ! for a single bin. Applies bam_sulfate_scale only to SULFATE (not volcanic).
  ! b4b operation order: (mmr * rho) first, then * ntm [* 2.0 for sulfate].
  !------------------------------------------------------------------------------
  subroutine get_bulk_num_and_mass(self, bin_ndx, ncol, rho, naer2, maerosol)
    class(bulk_aerosol_state), intent(in) :: self
    integer, intent(in)   :: bin_ndx
    integer, intent(in)   :: ncol
    real(r8), intent(in)  :: rho(:,:)       ! air density (kg/m3), (ncol,pver)
    real(r8), intent(out) :: naer2(:,:)     ! number concentration (#/m3)
    real(r8), intent(out) :: maerosol(:,:)  ! mass concentration (kg/m3)

    real(r8), pointer :: mmr(:,:)
    real(r8)          :: ntm
    character(len=32) :: aname

    call self%get_ambient_mmr(species_ndx=1, bin_ndx=bin_ndx, mmr=mmr)
    call rad_aer_get_props(self%list_idx_, bin_ndx, num_to_mass_aer=ntm, aername=aname)

    ! b4b operation order: (mmr * rho) first, then * ntm [* 2.0 for sulfate]
    !
    ! Note: only SULFATE gets the scale factor here.
    ! Volcanic aerosol (which also has spectype 'sulfate') does not get scaled in the
    ! ndrop_bam/CCN path (which only scales idxsul, SULFATE here)
    ! Ice nucleation has been unified to also use this path, but it does scale volcanic
    ! aerosol; it will apply this scale factor separately.
    maerosol(:ncol,:) = mmr(:ncol,:) * rho(:ncol,:)

    select case ( to_lower( aname(:4) ) )
    case ('sulf')
       naer2(:ncol,:) = maerosol(:ncol,:) * ntm * bam_sulfate_scale
    case default
       naer2(:ncol,:) = maerosol(:ncol,:) * ntm
    end select

  end subroutine get_bulk_num_and_mass

  ! NOTE on bit-for-bit: The base-class nuclice_get_numdens computes:
  !   size_wght * type_wght * num_col(#/kg) * rho * per_cm3
  ! where for BAM: num_col = mmr * ntm [* bam_sulfate_scale], size_wght = 1/25, type_wght = 1.0
  ! giving: (1/25) * 1.0 * (mmr * ntm) * rho * 1e-6
  !
  ! The original inline BAM code (nucleate_ice_cam.F90, removed) computed:
  !   naer2 = aer_mmr * rho * ntm  (mmr * rho first, then * ntm)
  !   dust_num = naer2 / 25 * 1e-6
  ! giving: (mmr * rho * ntm) / 25 * 1e-6
  !
  ! These differ only in floating-point operation order (associativity).
  ! It has been shown that this rearranging causes answer differences, so we
  ! use this subroutine to replicate the original behavior.
  subroutine nuclice_get_numdens_bam(self, aero_props, use_preexisting_ice, &
       ncol, nlev, rho, dust_num_col, sulf_num_col, soot_num_col, sulf_num_tot_col)
    !REMOVECAM: host-model specific dimensions
    use ppgrid, only: pcols, pver
    !REMOVECAM_END

    class(bulk_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    logical, intent(in) :: use_preexisting_ice
    integer, intent(in) :: ncol
    integer, intent(in) :: nlev
    real(r8), intent(in) :: rho(:,:)
    real(r8), intent(out) :: dust_num_col(:,:)
    real(r8), intent(out) :: sulf_num_col(:,:)
    real(r8), intent(out) :: soot_num_col(:,:)
    real(r8), intent(out) :: sulf_num_tot_col(:,:)

    real(r8) :: naer2_1bin(ncol,nlev)
    real(r8) :: maerosol_1bin(ncol,nlev)
    character(len=32) :: spectype, aname
    integer :: m, i, k
    real(r8), parameter :: per_cm3 = 1.e-6_r8

    dust_num_col(:,:) = 0._r8
    sulf_num_col(:,:) = 0._r8
    soot_num_col(:,:) = 0._r8
    sulf_num_tot_col(:,:) = 0._r8

    do m = 1, aero_props%nbins()
       call aero_props%species_type(m, 1, spectype)
       call self%get_bulk_num_and_mass(m, ncol, rho, naer2_1bin, maerosol_1bin)

       ! get_bulk_num_and_mass only applied bam_sulfate_scale to SULFATE (by name).
       ! For the nucleate_ice path, volcanic aerosol (spectype 'sulfate', name 'volc*')
       ! also needs the scale, matching the original inline code which scaled ALL
       ! spectype=='sulfate' bins including volcanic aerosol, so we will do it here:
       ! (but do not do it again for SULFATE)
       if (spectype == 'sulfate') then
          call rad_aer_get_props(self%list_idx_, m, aername=aname)
          if (to_lower(aname(:4)) == 'volc') then
             naer2_1bin(:ncol,:nlev) = naer2_1bin(:ncol,:nlev) * bam_sulfate_scale
          end if
       end if

       do k = 1, nlev
          do i = 1, ncol
             select case (trim(spectype))
             case ('dust')
                dust_num_col(i,k) = dust_num_col(i,k) + naer2_1bin(i,k) / 25._r8 * per_cm3
             case ('sulfate')
                sulf_num_col(i,k) = sulf_num_col(i,k) + naer2_1bin(i,k) / 25._r8 * per_cm3
                sulf_num_tot_col(i,k) = sulf_num_tot_col(i,k) + naer2_1bin(i,k) / 25._r8 * per_cm3
             case ('black-c')
                soot_num_col(i,k) = soot_num_col(i,k) + naer2_1bin(i,k) / 25._r8 * per_cm3
             end select
          end do
       end do
    end do

  end subroutine nuclice_get_numdens_bam

end module bulk_aerosol_state_mod
