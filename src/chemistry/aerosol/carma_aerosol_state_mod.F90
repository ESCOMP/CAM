module carma_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_state_mod, only: aerosol_state, ptr2d_t

  use rad_constituents, only: rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_num !, rad_cnst_get_bin_mmr
  use rad_constituents, only: rad_cnst_get_info_by_bin
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len

  use physconst, only: pi
  use carma_intr, only: carma_get_total_mmr, carma_get_dry_radius, carma_get_number, carma_get_number_cld
  use carma_intr, only: carma_get_group_by_name, carma_get_kappa, carma_get_dry_radius, carma_get_wet_radius
  use carma_intr, only: carma_get_wght_pct
  use ppgrid, only: begchunk, endchunk, pcols, pver

  implicit none

  private

  public :: carma_aerosol_state

  type, extends(aerosol_state) :: carma_aerosol_state
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
     procedure :: update_bin
     procedure :: hetfrz_size_wght
     procedure :: hygroscopicity
     procedure :: water_uptake
     procedure :: wgtpct
     procedure :: dry_volume
     procedure :: wet_volume
     procedure :: water_volume
     procedure :: wet_diameter

     final :: destructor

  end type carma_aerosol_state

  interface carma_aerosol_state
     procedure :: constructor
  end interface carma_aerosol_state

  real(r8), parameter :: four_thirds_pi = pi * 4._r8 / 3._r8

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf) result(newobj)
    type(physics_state), target, optional :: state
    type(physics_buffer_desc), pointer, optional :: pbuf(:)

    type(carma_aerosol_state), pointer :: newobj

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
    type(carma_aerosol_state), intent(inout) :: self

    nullify(self%state)
    nullify(self%pbuf)

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! sets transported components
  ! This aerosol model with the state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(carma_aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  ! This returns to current state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(carma_aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine get_transported

  !------------------------------------------------------------------------
  ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
  !------------------------------------------------------------------------
  function ambient_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    integer, intent(in) :: bin_ndx      ! bin index
    integer, intent(in) :: col_ndx      ! column index
    integer, intent(in) :: lyr_ndx      ! vertical layer index

    real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species

    real(r8) :: totmmr(pcols,pver)
    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call  carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call carma_get_total_mmr(self%state, igroup, ibin, totmmr, rc)

    mmr_tot = totmmr(col_ndx,lyr_ndx)

  end function ambient_total_bin_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_0list(self, species_ndx, bin_ndx, mmr)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    call rad_cnst_get_bin_mmr_by_idx(0, bin_ndx, species_ndx, 'a', self%state, self%pbuf, mmr)

  end subroutine get_ambient_mmr_0list

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given radiation diagnostics
  ! list index, species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_rlist(self, list_ndx, species_ndx, bin_ndx, mmr)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: list_ndx     ! rad climate list index
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    call rad_cnst_get_bin_mmr_by_idx(list_ndx, bin_ndx, species_ndx, 'a', self%state, self%pbuf, mmr)

  end subroutine get_ambient_mmr_rlist

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    call rad_cnst_get_bin_mmr_by_idx(0, bin_ndx, species_ndx, 'c', self%state, self%pbuf, mmr)

  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer :: num(:,:)      ! number mixing ratios

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr, ncol
    real(r8) :: nmr(pcols,pver)

    ncol = self%state%ncol

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call rad_cnst_get_bin_num(0, bin_ndx, 'a', self%state, self%pbuf, num)

    call carma_get_number(self%state, igroup, ibin, nmr, rc)

    num(:ncol,:) = nmr(:ncol,:)

  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer :: num(:,:)      ! number mixing ratios

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr, ncol
    real(r8) :: nmr(pcols,pver)

    ncol = self%state%ncol

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call  carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call rad_cnst_get_bin_num(0, bin_ndx, 'c', self%state, self%pbuf, num)

    call carma_get_number_cld(self%pbuf, igroup, ibin,  ncol, pver, nmr, rc)

    num(:ncol,:) = nmr(:ncol,:)

  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  ! returns interstitial and cloud-borne aerosol states
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(carma_aerosol_state), intent(in) :: self
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
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
    real(r8), intent(out) :: wght(:,:)

    character(len=aero_name_len) :: bin_name, shortname
    real(r8) :: rdry(ncol,nlev), rhopdry(ncol,nlev)
    integer :: i,k
    real(r8) :: diamdry
    integer :: igroup, ibin, rc, nchr

    wght = 0._r8

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)
    call  carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call carma_get_dry_radius(self%state, igroup, ibin, rdry, rhopdry, rc) ! m, kg/m3

    do k = 1,nlev
       do i = 1,ncol
          diamdry = rdry(i,k) * 2._r8 * 1.e6_r8 ! diameter in microns (from radius in m)
          if (diamdry >= 0.1_r8) then ! size threashold
             wght(i,k) = 1._r8
          end if
       end do
    end do

  end subroutine icenuc_size_wght_arr

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin, column and vertical layer
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_val(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
    real(r8), intent(out) :: wght

    real(r8) :: wght_arr(pcols,pver)

    call self%icenuc_size_wght(bin_ndx, self%state%ncol, pver, species_type, use_preexisting_ice, wght_arr)

    wght = wght_arr(col_ndx,lyr_ndx)

  end subroutine icenuc_size_wght_val

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine update_bin( self, bin_ndx, col_ndx, lyr_ndx, delmmr_sum, delnum_sum, tnd_ndx, dtime, tend )
    class(carma_aerosol_state), intent(in) :: self
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

    ! for updating num (num tendancies)
    ! -- nothing to do here for CARMA since num is calculated when needed

  end subroutine update_bin

  !------------------------------------------------------------------------------
  ! returns the volume-weighted fractions of aerosol subset `bin_ndx` that can act
  ! as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_size_wght(self, bin_ndx, ncol, nlev) result(wght)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels

    real(r8) :: wght(ncol,nlev)

    character(len=aero_name_len) :: bin_name, shortname
    real(r8) :: rdry(ncol,nlev), rhopdry(ncol,nlev)
    integer :: i,k
    real(r8) :: diamdry
    integer :: igroup, ibin, rc, nchr

    wght = 0._r8

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)
    call  carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call carma_get_dry_radius(self%state, igroup, ibin, rdry, rhopdry, rc) ! m, kg/m3

    do k = 1,nlev
       do i = 1,ncol
          diamdry = rdry(i,k) * 2._r8 * 1.e6_r8 ! diameter in microns (from radius in m)
          if (diamdry >= 0.1_r8) then ! size threashold
             wght(i,k) = 1._r8
          end if
       end do
    end do

  end function hetfrz_size_wght

  !------------------------------------------------------------------------------
  ! returns hygroscopicity for a given radiation diagnostic list number and
  ! bin number
  !------------------------------------------------------------------------------
  subroutine hygroscopicity(self, list_ndx, bin_ndx, kappa)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: list_ndx        ! rad climate list number
    integer, intent(in) :: bin_ndx         ! bin number
    real(r8), intent(out) :: kappa(:,:)    ! hygroscopicity (ncol,nlev)

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr, ncol

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call carma_get_kappa(self%state, igroup, ibin, kappa, rc)

  end subroutine hygroscopicity

  !------------------------------------------------------------------------------
  ! returns aerosol wet diameter and aerosol water concentration for a given
  ! radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  subroutine water_uptake(self, aero_props, list_idx, bin_idx, ncol, nlev, dgnumwet, qaerwat)

    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: list_idx             ! rad climate/diags list number
    integer, intent(in) :: bin_idx              ! bin number
    integer, intent(in) :: ncol                 ! number of columns
    integer, intent(in) :: nlev                 ! number of levels
    real(r8),intent(out) :: dgnumwet(ncol,nlev) ! aerosol wet diameter (m)
    real(r8),intent(out) :: qaerwat(ncol,nlev)  ! aerosol water concentration (g/g)

    dgnumwet = -huge(1._r8)
    qaerwat = -huge(1._r8)

  end subroutine water_uptake

  !------------------------------------------------------------------------------
  ! aerosol weight precent of H2SO4/H2O solution
  !------------------------------------------------------------------------------
  function wgtpct(self, ncol, nlev) result(wtp)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) ::  ncol, nlev
    real(r8) :: wtp(ncol,nlev)  ! weight precent of H2SO4/H2O solution for given icol, ilev

    wtp(:,:) = carma_get_wght_pct(ncol,nlev,self%state)

  end function wgtpct

  !------------------------------------------------------------------------------
  ! aerosol dry volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function dry_volume(self, aero_props, list_idx, bin_idx, ncol, nlev) result(vol)

    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_idx  ! rad climate/diags list number
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: raddry(pcols,pver)   ! dry radius (m)
    real(r8) :: rhodry(pcols,pver)   ! dry density (kg/m3)
    real(r8) :: nmr(pcols,pver)      ! number mixing ratio (#/kg)

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr

    call rad_cnst_get_info_by_bin(0, bin_idx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call  carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    vol = 0._r8

    call carma_get_dry_radius(self%state, igroup, ibin, raddry, rhodry, rc)
    call carma_get_number(self%state, igroup, ibin, nmr, rc)

    vol(:ncol,:) = four_thirds_pi * (raddry(:ncol,:)**3) * nmr(:ncol,:) ! units = m3/kg

  end function dry_volume

  !------------------------------------------------------------------------------
  ! aerosol wet volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function wet_volume(self, aero_props, list_idx, bin_idx, ncol, nlev) result(vol)

    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_idx  ! rad climate/diags list number
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: radwet(pcols,pver)   ! wet radius (m)
    real(r8) :: rhowet(pcols,pver)   ! wet density (kg/m3)
    real(r8) :: nmr(pcols,pver)      ! number mixing ratio (#/kg)

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr

    call rad_cnst_get_info_by_bin(0, bin_idx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call  carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    vol = 0._r8

    call carma_get_wet_radius(self%state, igroup, ibin, radwet, rhowet, rc)
    call carma_get_number(self%state, igroup, ibin, nmr, rc)

    vol(:ncol,:) = four_thirds_pi * (radwet(:ncol,:)**3) * nmr(:ncol,:) ! units = m3/kg

  end function wet_volume

  !------------------------------------------------------------------------------
  ! aerosol water volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function water_volume(self, aero_props, list_idx, bin_idx, ncol, nlev) result(vol)

    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_idx  ! rad climate/diags list number
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: wetvol(ncol,nlev)
    real(r8) :: dryvol(ncol,nlev)

    wetvol = self%wet_volume(aero_props, list_idx, bin_idx, ncol, nlev)
    dryvol = self%dry_volume(aero_props, list_idx, bin_idx, ncol, nlev)

    vol(:ncol,:) = wetvol(:ncol,:) - dryvol(:ncol,:)

    where (vol<0._r8)
       vol = 0._r8
    end where

  end function water_volume

  !------------------------------------------------------------------------------
  ! aerosol wet diameter
  !------------------------------------------------------------------------------
  function wet_diameter(self, bin_idx, ncol, nlev) result(diam)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_idx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: diam(ncol,nlev)

    real(r8) :: radwet(pcols,pver)   !! wet radius (m)
    real(r8) :: rhowet(pcols,pver)   !! wet density (kg/m3)

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr

    call rad_cnst_get_info_by_bin(0, bin_idx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call  carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call carma_get_wet_radius(self%state, igroup, ibin, radwet, rhowet, rc)

    diam(:ncol,:nlev) = 2._r8*radwet(:ncol,:nlev)

  end function wet_diameter

end module carma_aerosol_state_mod
