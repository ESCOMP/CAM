!===============================================================================
! Dust emissions for the Modal Aerosol Model
!   Portable science routines split from modal_aero/dust_model.F90: selection
!   of the emitted dust size distribution, and rebinning of the coupler dust
!   flux into modal mass and number surface fluxes.
!   Host constants and index maps are passed as arguments; array sizing is by
!   runtime ncol.
!===============================================================================
module modal_dust_emissions
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: modal_dust_emissions_init
  public :: modal_dust_emissions_run

contains

  !=============================================================================
  ! Set the emitted dust size distribution (per-bin mass fractions) and the
  ! bin mass-weighted diameters used for the number flux conversion.
  !=============================================================================
  subroutine modal_dust_emissions_init( ntot_amode, dust_nbin, pi, rair, gravit, &
                                        dust_emis_sclfctr, dust_dmt_vwr, errmsg, errflg )
    use dust_common, only: dust_set_params

    integer,  intent(in)  :: ntot_amode            ! number of aerosol modes
    integer,  intent(in)  :: dust_nbin             ! number of dust bins (mass species)
    real(r8), intent(in)  :: pi                    ! host model constants
    real(r8), intent(in)  :: rair                  ! gas constant for dry air (J/K/kg)
    real(r8), intent(in)  :: gravit                ! gravitational acceleration (m/s2)
    real(r8), intent(out) :: dust_emis_sclfctr(:)  ! mass fraction of emissions per bin
    real(r8), intent(out) :: dust_dmt_vwr(:)       ! mass-weighted diameter per bin (m)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(r8) :: dust_dmt_grd(dust_nbin+1)  ! bin diameter edges (m)
    real(r8) :: dust_stk_crc(dust_nbin)    ! Stokes correction from dust_set_params; unused by emissions

    errmsg = ''
    errflg = 0

    ! dmleung edited the mass fraction of the emitted dust size distribution. 27 Oct 2025 ++
    ! The new mass fraction comes from Jun Meng et al. (2022) and MERRA-2.
    ! Jun Meng's table indicates 2.1 % mass for 0.1-1 um and 97.9 % mass for 1-10 um.
    ! ref: https://zenodo.org/records/6344524
    ! MERRA-2 dust emissions indicate 6 % mass for 0.1-1 um (bin1) and 94 % for 1-10 um (bin2-5).
    ! dmleung adopts 2.1 % mass for 0.1-1 um and 97.9 % mass for 1-10 um for dust.
    ! Distributing more mass to accumulation mode allows a longer lifetime of dust, reducing
    ! low dust biases over remote oceans and reducing high dust biases over the Sahara.
    ! This change impacts both Zender_2003 dust and Leung_2023 dust.
    if ( ntot_amode == 3 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.021_r8,0.979_r8 /)
    elseif ( ntot_amode == 4 .or. ntot_amode == 5 ) then
       dust_dmt_grd(:) = (/ 0.01e-6_r8, 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /)
       dust_emis_sclfctr(:) = (/ 1.65E-05_r8, 0.021_r8, 0.979_r8 /)
    else if( ntot_amode == 7 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.12_r8, 0.88_r8 /)
    endif
    ! dmleung --

    call dust_set_params( nbin=dust_nbin, dmt_grd=dust_dmt_grd,      &
                          dmt_vwr=dust_dmt_vwr, stk_crc=dust_stk_crc, &
                          pi=pi, rair=rair, gravit=gravit,            &
                          errmsg=errmsg, errflg=errflg )

  end subroutine modal_dust_emissions_init

  !===============================================================================
  ! Rebin and adjust the incoming coupler dust flux into per-bin mass and
  ! number surface fluxes.
  !===============================================================================
  subroutine modal_dust_emissions_run( ncol, dust_nbin, dust_indices, dust_emis_sclfctr, &
                                       dust_dmt_vwr, dust_emis_fact,                     &
                                       zender_soil_erod_from_atm, soil_erodibility,      &
                                       dust_flux_in, pi, cflx, soil_erod )
    use dust_common, only: dust_density

  ! args
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: dust_nbin             ! number of dust bins (mass species)
    integer,  intent(in)    :: dust_indices(:)       ! constituent indices: mass bins, then number bins
    real(r8), intent(in)    :: dust_emis_sclfctr(:)  ! mass fraction of emissions per bin
    real(r8), intent(in)    :: dust_dmt_vwr(:)       ! mass-weighted diameter per bin (m)
    real(r8), intent(in)    :: dust_emis_fact        ! tuning parameter for dust emissions
    logical,  intent(in)    :: zender_soil_erod_from_atm ! Zender_2003 with soil erodibility applied in atm
    real(r8), intent(in)    :: soil_erodibility(:)   ! soil erodibility factor (used only when
                                                     ! zender_soil_erod_from_atm)
    real(r8), intent(in)    :: dust_flux_in(:,:)     ! dust fluxes from the coupler (kg/m2/s, negative down)
    real(r8), intent(in)    :: pi
    real(r8), intent(inout) :: cflx(:,:)             ! constituent surface fluxes (kg/m2/s)
    real(r8), intent(out)   :: soil_erod(:)          ! thresholded soil erodibility
                                                     ! (not set on the Leung_2023 branch)

  ! local vars
    integer :: i, m, idst, inum
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8

    ! set dust emissions

    if (zender_soil_erod_from_atm) then   ! Zender_2003 dust emissions
       col_loop1: do i = 1,ncol
          soil_erod(i) = soil_erodibility(i)
          if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

          ! rebin and adjust dust emissons.
          do m = 1,dust_nbin
             idst = dust_indices(m)
             cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
                  * dust_emis_sclfctr(m)*soil_erod(i)/dust_emis_fact*1.15_r8
             x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))
             inum = dust_indices(m+dust_nbin)
             cflx(i,inum) = cflx(i,idst)*x_mton
          enddo
       enddo col_loop1
    else ! Leung_2023 dust emissions

       col_loop2: do i = 1,ncol
          ! rebin and adjust dust emissons.
          do m = 1,dust_nbin
             idst = dust_indices(m)

             cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
                  * dust_emis_sclfctr(m) / dust_emis_fact
             x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))
             inum = dust_indices(m+dust_nbin)
             cflx(i,inum) = cflx(i,idst)*x_mton
          enddo
       enddo col_loop2
    end if

  end subroutine modal_dust_emissions_run

end module modal_dust_emissions
