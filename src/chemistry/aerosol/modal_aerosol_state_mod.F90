module modal_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_spfn_mod, only: erf => shr_spfn_erf
  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use rad_constituents, only: rad_cnst_get_aer_mmr, rad_cnst_get_mode_num, rad_cnst_get_info
  use rad_constituents, only: rad_cnst_get_mode_props
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use physconst,  only: rhoh2o

  implicit none

  private

  public :: modal_aerosol_state

  type, extends(aerosol_state) :: modal_aerosol_state
     private
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
   contains

     procedure :: get_transported
     procedure :: set_transported
     procedure :: ambient_total_bin_mmr
     procedure :: get_ambient_mmr_0list
     procedure :: get_ambient_mmr_rlist
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

     final :: destructor

  end type modal_aerosol_state

  interface modal_aerosol_state
     procedure :: constructor
  end interface modal_aerosol_state

  real(r8), parameter :: rh2odens = 1._r8/rhoh2o

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf) result(newobj)
    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(modal_aerosol_state), pointer :: newobj

    integer :: ierr

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
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
  ! sets transported components
  ! This aerosol model with the state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(modal_aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  ! This returns to current state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(modal_aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine get_transported

  !------------------------------------------------------------------------
  ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
  !------------------------------------------------------------------------
  function ambient_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: bin_ndx      ! bin index
    integer, intent(in) :: col_ndx      ! column index
    integer, intent(in) :: lyr_ndx      ! vertical layer index

    real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species
    real(r8),pointer :: mmrptr(:,:)
    integer :: spec_ndx

    mmr_tot = 0._r8

    do spec_ndx=1,aero_props%nspecies(bin_ndx)
       call rad_cnst_get_aer_mmr(0, bin_ndx, spec_ndx, 'a', self%state, self%pbuf, mmrptr)
       mmr_tot = mmr_tot + mmrptr(col_ndx,lyr_ndx)
    end do

  end function ambient_total_bin_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_0list(self, species_ndx, bin_ndx, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    call rad_cnst_get_aer_mmr(0, bin_ndx, species_ndx, 'a', self%state, self%pbuf, mmr)
  end subroutine get_ambient_mmr_0list

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given radiation diagnostics
  ! list index, species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_rlist(self, list_ndx, species_ndx, bin_ndx, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: list_ndx     ! rad climate list index
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    call rad_cnst_get_aer_mmr(list_ndx, bin_ndx, species_ndx, 'a', self%state, self%pbuf, mmr)
  end subroutine get_ambient_mmr_rlist

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

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

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_arr(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
    real(r8), intent(out) :: wght(:,:)

    character(len=aero_name_len) :: modetype
    real(r8), pointer :: dgnum(:,:,:)    ! mode dry radius
    real(r8) :: sigmag_aitken
    integer :: i,k

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)

    wght = 0._r8

    select case ( trim(species_type) )
    case('dust')
       if (modetype=='coarse' .or. modetype=='coarse_dust') then
          wght(:ncol,:) = 1._r8
       end if
    case('sulfate')
       if (modetype=='aitken') then
          if ( use_preexisting_ice ) then
             wght(:ncol,:) = 1._r8
          else
             call rad_cnst_get_mode_props(0, bin_ndx, sigmag=sigmag_aitken)
             call pbuf_get_field(self%pbuf, pbuf_get_index('DGNUM' ), dgnum)
             do k = 1,nlev
                do i = 1,ncol
                   if (dgnum(i,k,bin_ndx) > 0._r8) then
                      ! only allow so4 with D>0.1 um in ice nucleation
                      wght(i,k) = max(0._r8,(0.5_r8 - 0.5_r8* &
                           erf(log(0.1e-6_r8/dgnum(i,k,bin_ndx))/ &
                           (2._r8**0.5_r8*log(sigmag_aitken)))  ))
                   end if
                end do
             end do
          endif
       endif
    case('black-c')
       if (modetype=='accum') then
          wght(:ncol,:) = 1._r8
       endif
    case('sulfate_strat')
       if (modetype=='accum' .or. modetype=='coarse' .or. modetype=='coarse_strat') then
          wght(:ncol,:) = 1._r8
       endif
    end select

  end subroutine icenuc_size_wght_arr

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin, column and vertical layer
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_val(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
    real(r8), intent(out) :: wght

    character(len=aero_name_len) :: modetype
    real(r8), pointer :: dgnum(:,:,:)    ! mode dry radius
    real(r8) :: sigmag_aitken

    wght = 0._r8

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)

    select case ( trim(species_type) )
    case('dust')
       if (modetype=='coarse' .or. modetype=='coarse_dust') then
          wght = 1._r8
       end if
    case('sulfate')
       if (modetype=='aitken') then
          if ( use_preexisting_ice ) then
             wght = 1._r8
          else
             call rad_cnst_get_mode_props(0, bin_ndx, sigmag=sigmag_aitken)
             call pbuf_get_field(self%pbuf, pbuf_get_index('DGNUM' ), dgnum)

             if (dgnum(col_ndx,lyr_ndx,bin_ndx) > 0._r8) then
                ! only allow so4 with D>0.1 um in ice nucleation
                wght = max(0._r8,(0.5_r8 - 0.5_r8* &
                     erf(log(0.1e-6_r8/dgnum(col_ndx,lyr_ndx,bin_ndx))/ &
                     (2._r8**0.5_r8*log(sigmag_aitken)))  ))

             end if
          endif
       endif
    case('black-c')
       if (modetype=='accum') then
          wght = 1._r8
       endif
    case('sulfate_strat')
       if (modetype=='accum' .or. modetype=='coarse' .or. modetype=='coarse_strat') then
          wght = 1._r8
       endif
    end select

  end subroutine icenuc_size_wght_val

  !------------------------------------------------------------------------------
  ! returns aerosol type weights for a given aerosol type and bin
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, wght, cloud_borne)

    use aerosol_properties_mod, only: aerosol_properties

    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wght(:,:)            ! type weights
    logical, optional, intent(in) :: cloud_borne  ! if TRUE cloud-borne aerosols are used
                                                  ! otherwise ambient aerosols are used

    character(len=aero_name_len) :: modetype

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)

    wght = 0._r8

    if (species_type == 'dust') then
       if (modetype=='coarse_dust') then
          wght(:ncol,:) = 1._r8
       else
          call self%icenuc_type_wght_base(bin_ndx, ncol, nlev, species_type, aero_props, rho, wght, cloud_borne)
       end if
    else if (species_type == 'sulfate_strat') then
       if (modetype=='accum') then
          wght(:ncol,:) = 1._r8
       elseif ( modetype=='coarse' .or. modetype=='coarse_strat') then
          call self%icenuc_type_wght_base(bin_ndx, ncol, nlev, species_type, aero_props, rho, wght, cloud_borne)
       endif
    else
       wght(:ncol,:) = 1._r8
    end if

  end subroutine icenuc_type_wght

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine update_bin( self, bin_ndx, col_ndx, lyr_ndx, delmmr_sum, delnum_sum, tnd_ndx, dtime, tend )
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    real(r8),intent(in) :: delmmr_sum             ! mass mixing ratio change summed over all species in bin
    real(r8),intent(in) :: delnum_sum             ! number mixing ratio change summed over all species in bin
    integer, intent(in) :: tnd_ndx                ! tendency index
    real(r8),intent(in) :: dtime                  ! time step size (sec)
    real(r8),intent(inout) :: tend(:,:,:)         ! tendency

    real(r8), pointer :: amb_num(:,:)
    real(r8), pointer :: cld_num(:,:)

    call self%get_ambient_num(bin_ndx, amb_num)
    call self%get_cldbrne_num(bin_ndx, cld_num)

    ! if there is no bin mass compute updates/tendencies for bin number
    ! -- apply the total number change to bin number
    if (tnd_ndx>0) then
       tend(col_ndx,lyr_ndx,tnd_ndx) = -delnum_sum/dtime
    else
       amb_num(col_ndx,lyr_ndx) = amb_num(col_ndx,lyr_ndx) - delnum_sum
    end if

    ! apply the total number change to bin number
    cld_num(col_ndx,lyr_ndx) = cld_num(col_ndx,lyr_ndx) + delnum_sum

  end subroutine update_bin

  !------------------------------------------------------------------------------
  ! returns the volume-weighted fractions of aerosol subset `bin_ndx` that can act
  ! as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_size_wght(self, bin_ndx, ncol, nlev) result(wght)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels

    real(r8) :: wght(ncol,nlev)

    character(len=aero_name_len) :: modetype

    wght(:,:) = 1._r8

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)

    if (trim(modetype) == 'aitken') then
       wght(:,:) = 0._r8
    end if

  end function hetfrz_size_wght

  !------------------------------------------------------------------------------
  ! returns hygroscopicity for a given radiation diagnostic list number and
  ! bin number
  !------------------------------------------------------------------------------
  function hygroscopicity(self, list_ndx, bin_ndx) result(kappa)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: list_ndx        ! rad climate list number
    integer, intent(in) :: bin_ndx         ! bin number

    real(r8), pointer :: kappa(:,:)        ! hygroscopicity (ncol,nlev)

    nullify(kappa)

  end function hygroscopicity

  !------------------------------------------------------------------------------
  ! returns aerosol wet diameter and aerosol water concentration for a given
  ! radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  subroutine water_uptake(self, aero_props, list_idx, bin_idx, ncol, nlev, dgnumwet, qaerwat)
    use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
    use modal_aero_calcsize,    only: modal_aero_calcsize_diag

    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: list_idx             ! rad climate/diags list number
    integer, intent(in) :: bin_idx              ! bin number
    integer, intent(in) :: ncol                 ! number of columns
    integer, intent(in) :: nlev                 ! number of levels
    real(r8),intent(out) :: dgnumwet(ncol,nlev) ! aerosol wet diameter (m)
    real(r8),intent(out) :: qaerwat(ncol,nlev)  ! aerosol water concentration (g/g)

    integer :: istat, nmodes
    real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
    real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
    real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
    real(r8), pointer :: wetdens_m(:,:,:)  !
    real(r8), pointer :: hygro_m(:,:,:)  !
    real(r8), pointer :: dryvol_m(:,:,:)  !
    real(r8), pointer :: dryrad_m(:,:,:)  !
    real(r8), pointer :: drymass_m(:,:,:)  !
    real(r8), pointer :: so4dryvol_m(:,:,:)  !
    real(r8), pointer :: naer_m(:,:,:)  !

    nmodes = aero_props%nbins()

    if (list_idx == 0) then
       ! water uptake and wet radius for the climate list has already been calculated
       call pbuf_get_field(self%pbuf, pbuf_get_index('DGNUMWET'), dgnumwet_m)
       call pbuf_get_field(self%pbuf, pbuf_get_index('QAERWAT'),  qaerwat_m)

       dgnumwet(:ncol,:nlev) = dgnumwet_m(:ncol,:nlev,bin_idx)
       qaerwat (:ncol,:nlev) =  qaerwat_m(:ncol,:nlev,bin_idx)

    else
       ! If doing a diagnostic calculation then need to calculate the wet radius
       ! and water uptake for the diagnostic modes
       allocate(dgnumdry_m(ncol,nlev,nmodes),  dgnumwet_m(ncol,nlev,nmodes), &
                qaerwat_m(ncol,nlev,nmodes),   wetdens_m(ncol,nlev,nmodes), &
                hygro_m(ncol,nlev,nmodes),     dryvol_m(ncol,nlev,nmodes), &
                dryrad_m(ncol,nlev,nmodes),    drymass_m(ncol,nlev,nmodes),  &
                so4dryvol_m(ncol,nlev,nmodes), naer_m(ncol,nlev,nmodes), stat=istat)
       if (istat > 0) then
          dgnumwet = -huge(1._r8)
          qaerwat = -huge(1._r8)
          return
       end if
       call modal_aero_calcsize_diag(self%state, self%pbuf, list_idx, dgnumdry_m, hygro_m, &
                                     dryvol_m, dryrad_m, drymass_m, so4dryvol_m, naer_m)
       call modal_aero_wateruptake_dr(self%state, self%pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                      qaerwat_m, wetdens_m,  hygro_m, dryvol_m, dryrad_m, &
                                      drymass_m, so4dryvol_m, naer_m)

       dgnumwet(:ncol,:nlev) = dgnumwet_m(:ncol,:nlev,bin_idx)
       qaerwat (:ncol,:nlev) =  qaerwat_m(:ncol,:nlev,bin_idx)

       deallocate(dgnumdry_m)
       deallocate(dgnumwet_m)
       deallocate(qaerwat_m)
       deallocate(wetdens_m)
       deallocate(hygro_m)
       deallocate(dryvol_m)
       deallocate(dryrad_m)
       deallocate(drymass_m)
       deallocate(so4dryvol_m)
       deallocate(naer_m)
    endif


  end subroutine water_uptake

  !------------------------------------------------------------------------------
  ! aerosol dry volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function dry_volume(self, aero_props, list_idx, bin_idx, ncol, nlev) result(vol)

    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_idx  ! rad climate/diags list number
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8), pointer :: mmr(:,:)
    real(r8) :: specdens              ! species density (kg/m3)

    integer :: ispec

    vol(:,:) = 0._r8

    do ispec = 1, aero_props%nspecies(list_idx,bin_idx)
       call self%get_ambient_mmr(list_idx, ispec, bin_idx, mmr)
       call aero_props%get(bin_idx, ispec, list_ndx=list_idx, density=specdens)
       vol(:ncol,:) = vol(:ncol,:) + mmr(:ncol,:)/specdens
    end do

  end function dry_volume

  !------------------------------------------------------------------------------
  ! aerosol wet volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function wet_volume(self, aero_props, list_idx, bin_idx, ncol, nlev) result(vol)

    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_idx  ! rad climate/diags list number
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: dryvol(ncol,nlev)
    real(r8) :: watervol(ncol,nlev)

    dryvol = self%dry_volume(aero_props, list_idx, bin_idx, ncol, nlev)
    watervol = self%water_volume(aero_props, list_idx, bin_idx, ncol, nlev)

    vol = watervol + dryvol

  end function wet_volume

  !------------------------------------------------------------------------------
  ! aerosol water volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function water_volume(self, aero_props, list_idx, bin_idx, ncol, nlev) result(vol)

    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_idx  ! rad climate/diags list number
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: dgnumwet(ncol,nlev)
    real(r8) :: qaerwat(ncol,nlev)

    call self%water_uptake(aero_props, list_idx, bin_idx, ncol, nlev, dgnumwet, qaerwat)

    vol(:ncol,:nlev) = qaerwat(:ncol,:nlev)*rh2odens
    where (vol<0._r8)
       vol = 0._r8
    end where

  end function water_volume

end module modal_aerosol_state_mod
