!===============================================================================
! Sea salt emissions for the Modal Aerosol Model
!   Portable science routine split from modal_aero/seasalt_model.F90, plus the
!   10 m wind extrapolation moved from modal_aero/aero_model.F90: sea salt
!   section number fluxes accumulated into modal number and mass surface
!   fluxes over the ocean fraction.
!   Host constants and index maps are passed as arguments; array sizing is by
!   runtime ncol.
!===============================================================================
module modal_seasalt_emissions
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: modal_seasalt_emissions_run

  ! Sea salt aerosol material density used for the number-to-mass flux
  ! conversion (value from CAM mo_constants).
  real(r8), parameter :: seasalt_density = 2.2e+3_r8  ! [kg m-3]

contains

  !=============================================================================
  !=============================================================================
  subroutine modal_seasalt_emissions_run( ncol, nslt, seasalt_indices, emis_scale, &
                                          u_bottom, v_bottom, zmid_bottom,         &
                                          srf_temp, ocnfrc, pi, cflx )

    use sslt_sections, only: nsections, fluxes, Dg, rdry

    ! dummy arguments
    integer,  intent(in) :: ncol
    integer,  intent(in) :: nslt                ! number of seasalt bins (mass species)
    integer,  intent(in) :: seasalt_indices(:)  ! constituent indices: mass bins, then number bins
    real(r8), intent(in) :: emis_scale          ! sea salt emission tuning factor
    real(r8), intent(in) :: u_bottom(:)         ! bottom layer zonal wind (m/s)
    real(r8), intent(in) :: v_bottom(:)         ! bottom layer meridional wind (m/s)
    real(r8), intent(in) :: zmid_bottom(:)      ! bottom layer midpoint geopotential height above surface (m)
    real(r8), intent(in) :: srf_temp(:)         ! sea surface temperature (K)
    real(r8), intent(in) :: ocnfrc(:)           ! ocean fraction
    real(r8), intent(in) :: pi
    real(r8), intent(inout) :: cflx(:,:)        ! constituent surface fluxes (kg/m2/s or #/m2/s)

    ! local vars
    integer  :: mn, mm, ibin, i
    real(r8) :: fi(ncol,nsections)
    real(r8) :: u10cubed(ncol)
    real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model

    real(r8) :: sst_sz_range_lo (nslt)
    real(r8) :: sst_sz_range_hi (nslt)

    u10cubed(:ncol)=sqrt(u_bottom(:ncol)**2+v_bottom(:ncol)**2)

    ! move the winds to 10m high from the midpoint of the gridbox:
    ! follows Tie and Seinfeld and Pandis, p.859 with math.

    u10cubed(:ncol)=u10cubed(:ncol)*log(10._r8/z0)/log(zmid_bottom(:ncol)/z0)

    ! we need them to the 3.41 power, according to Gong et al., 1997:
    u10cubed(:ncol)=u10cubed(:ncol)**3.41_r8

    if (nslt==4) then
       sst_sz_range_lo (:) = (/ 0.08e-6_r8, 0.02e-6_r8, 0.3e-6_r8,  1.0e-6_r8 /) ! accu, aitken, fine, coarse
       sst_sz_range_hi (:) = (/ 0.3e-6_r8,  0.08e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /)
    else if (nslt==3) then
       sst_sz_range_lo (:) =  (/ 0.08e-6_r8,  0.02e-6_r8,  1.0e-6_r8 /)  ! accu, aitken, coarse
       sst_sz_range_hi (:) =  (/ 1.0e-6_r8,   0.08e-6_r8, 10.0e-6_r8 /)
    endif

    fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

    do ibin = 1,nslt
       mm = seasalt_indices(ibin)
       mn = seasalt_indices(nslt+ibin)

       if (mn>0) then
          do i=1, nsections
             if (Dg(i).ge.sst_sz_range_lo(ibin) .and. Dg(i).lt.sst_sz_range_hi(ibin)) then
                cflx(:ncol,mn)=cflx(:ncol,mn)+fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  !++ ag: scale sea-salt
             endif
          enddo
       endif

       cflx(:ncol,mm)=0.0_r8
       do i=1, nsections
          if (Dg(i).ge.sst_sz_range_lo(ibin) .and. Dg(i).lt.sst_sz_range_hi(ibin)) then
             cflx(:ncol,mm)=cflx(:ncol,mm)+fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  &   !++ ag: scale sea-salt
                  *4._r8/3._r8*pi*rdry(i)**3*seasalt_density  ! should use dry size, convert from number to mass flux (kg/m2/s)
          endif
       enddo

    enddo

  end subroutine modal_seasalt_emissions_run

end module modal_seasalt_emissions
