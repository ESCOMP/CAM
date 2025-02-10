!-------------------------------------------------------------------------------
! This module uses the solar irradiance data 
! to provide a spectral scaling factor
! to approximate the spectral distribution of irradiance
! when the radiation scheme might use a different solar source function
!-------------------------------------------------------------------------------
module rad_solar_var

  use shr_kind_mod ,     only : r8 => shr_kind_r8
  use radconstants,      only : nswbands, get_sw_spectral_boundaries, band2gpt_sw
  use solar_irrad_data,  only : sol_irrad, we, nbins, has_spectrum, sol_tsi
  use solar_irrad_data,  only : do_spctrl_scaling
  use cam_abortutils,    only : endrun
  use error_messages,    only : alloc_err

  implicit none
  save

  private
  public :: rad_solar_var_init
  public :: get_variability

  real(r8), allocatable :: irrad(:)           ! solar irradiance at model timestep in each band

  real(r8), allocatable :: radbinmax(:)
  real(r8), allocatable :: radbinmin(:)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

  subroutine rad_solar_var_init( )

    integer :: ierr
    integer :: radmax_loc

    if ( do_spctrl_scaling ) then

       if ( .not.has_spectrum ) then
          call endrun('rad_solar_var_init: solar input file must have irradiance spectrum')
       endif

       allocate (radbinmax(nswbands),stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for radbinmax')
       end if

       allocate (radbinmin(nswbands),stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for radbinmin')
       end if

       allocate (irrad(nswbands), stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for irrad')
       end if

       call get_sw_spectral_boundaries(radbinmin, radbinmax, 'nm')

       ! Make sure that the far-IR is included, even if radiation grid does not
       ! extend that far down. 10^5 nm corresponds to a wavenumber of
       ! 100 cm^-1.
       radmax_loc = maxloc(radbinmax,1)
       radbinmax(radmax_loc) = max(100000._r8,radbinmax(radmax_loc))

    endif

  end subroutine rad_solar_var_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  subroutine get_variability(toa_flux, sfac) 

     ! Arguments 
     real(r8), intent(in)  :: toa_flux(:,:) ! TOA flux to be scaled (columns,gpts)
     real(r8), intent(out) :: sfac(:,:)     ! scaling factors (columns,gpts)

     ! Local variables 
     integer :: i, j, istat, gpt_start, gpt_end, ncols
     real(r8), allocatable :: scale(:)
     character(len=*), parameter :: sub = 'get_variability'
    
     if (do_spctrl_scaling) then 

        ! Determine target irradiance for each band
        call integrate_spectrum(nbins, nswbands, we, radbinmin, radbinmax, sol_irrad, irrad)

        ncols = size(toa_flux, 1)
        allocate(scale(ncols), stat=istat)
        call alloc_err(istat, sub, 'scale', ncols)

        do i = 1, nswbands 
           gpt_start = band2gpt_sw(1,i) 
           gpt_end   = band2gpt_sw(2,i) 
           scale = spread(irrad(i), 1, ncols) / sum(toa_flux(:, gpt_start:gpt_end), dim=2)
           do j = gpt_start, gpt_end
              sfac(:,j) = scale
           end do
        end do

     else 
        sfac(:,:) = sol_tsi / spread(sum(toa_flux, 2), 2, size(toa_flux, 2))
     end if
  end subroutine get_variability


!-------------------------------------------------------------------------------
! private method.........
!-------------------------------------------------------------------------------

  subroutine integrate_spectrum( nsrc, ntrg, src_x, min_trg, max_trg, src, trg )

    use mo_util, only : rebin

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)  :: nsrc                  ! dimension source array
    integer,  intent(in)  :: ntrg                  ! dimension target array
    real(r8), intent(in)  :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)  :: max_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: min_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: src(nsrc)             ! source array
    real(r8), intent(out) :: trg(ntrg)             ! target array
 
    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    real(r8) :: trg_x(2), targ(1)         ! target coordinates
    integer  :: i

    do i = 1, ntrg

       trg_x(1) = min_trg(i)
       trg_x(2) = max_trg(i)

       call rebin( nsrc, 1, src_x, trg_x, src(1:nsrc), targ(:) )
       ! W/m2/nm --> W/m2
       trg( i ) = targ(1)*(trg_x(2)-trg_x(1))

    enddo


  end subroutine integrate_spectrum

end module rad_solar_var
