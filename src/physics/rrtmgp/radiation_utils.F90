module radiation_utils
  use ccpp_kinds,       only: kind_phys
  use interpolate_data, only: interp_type, lininterp_init, lininterp, &
                              extrap_method_bndry

  public :: radiation_utils_init
  public :: get_sw_spectral_boundaries_ccpp
  public :: get_lw_spectral_boundaries_ccpp
  public :: get_mu_lambda_weights_ccpp

  real(kind_phys), allocatable :: wavenumber_low_shortwave(:)
  real(kind_phys), allocatable :: wavenumber_high_shortwave(:)
  real(kind_phys), allocatable :: wavenumber_low_longwave(:)
  real(kind_phys), allocatable :: wavenumber_high_longwave(:)
  integer :: nswbands
  integer :: nlwbands
  logical :: wavenumber_boundaries_set = .false.

contains

  subroutine radiation_utils_init(nswbands_in, nlwbands_in, low_shortwave, high_shortwave, &
                  low_longwave, high_longwave, errmsg, errflg)
    integer, intent(in) :: nswbands_in
    integer, intent(in) :: nlwbands_in
    real(kind_phys), intent(in) :: low_shortwave(:)
    real(kind_phys), intent(in) :: high_shortwave(:)
    real(kind_phys), intent(in) :: low_longwave(:)
    real(kind_phys), intent(in) :: high_longwave(:)
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg
    character(len=256) :: alloc_errmsg

    errflg = 0
    errmsg = ''
    nswbands = nswbands_in
    nlwbands = nlwbands_in
    allocate(wavenumber_low_shortwave(nswbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_low_shortwave, message: ', &
          alloc_errmsg
    end if
    allocate(wavenumber_high_shortwave(nswbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_high_shortwave, message: ', &
          alloc_errmsg
    end if
    allocate(wavenumber_low_longwave(nlwbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_low_longwave, message: ', &
          alloc_errmsg
    end if
    allocate(wavenumber_high_longwave(nlwbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg,'(a,a)') 'radiation_utils_init: failed to allocate wavenumber_high_longwave, message: ', &
          alloc_errmsg
    end if

    wavenumber_low_shortwave = low_shortwave
    wavenumber_high_shortwave = high_shortwave
    wavenumber_low_longwave = low_longwave
    wavenumber_high_longwave = high_longwave

    wavenumber_boundaries_set = .true.

  end subroutine radiation_utils_init

!=========================================================================================

 subroutine get_sw_spectral_boundaries_ccpp(low_boundaries, high_boundaries, units, errmsg, errflg)

   ! provide spectral boundaries of each shortwave band

   real(kind_phys), dimension(:), intent(out) :: low_boundaries
   real(kind_phys), dimension(:), intent(out) :: high_boundaries
   character(*), intent(in) :: units ! requested units
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   character(len=*), parameter :: sub = 'get_sw_spectral_boundaries_ccpp'
   !----------------------------------------------------------------------------

   ! Set error variables
   errmsg = ''
   errflg = 0

   if (.not. wavenumber_boundaries_set) then
      write(errmsg,'(a,a)') sub, ': ERROR, wavenumber boundaries not set.'
   end if

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries = wavenumber_low_shortwave
      high_boundaries = wavenumber_high_shortwave
   case('m','meter','meters')
      low_boundaries = 1.e-2_kind_phys/wavenumber_high_shortwave
      high_boundaries = 1.e-2_kind_phys/wavenumber_low_shortwave
   case('nm','nanometer','nanometers')
      low_boundaries = 1.e7_kind_phys/wavenumber_high_shortwave
      high_boundaries = 1.e7_kind_phys/wavenumber_low_shortwave
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries = 1.e4_kind_phys/wavenumber_high_shortwave
      high_boundaries = 1.e4_kind_phys/wavenumber_low_shortwave
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._kind_phys/wavenumber_high_shortwave
      high_boundaries = 1._kind_phys/wavenumber_low_shortwave
   case default
      write(errmsg, '(a,a,a)') sub, ': ERROR, requested spectral units not recognized: ', units
      errflg = 1
   end select

 end subroutine get_sw_spectral_boundaries_ccpp

!=========================================================================================

subroutine get_lw_spectral_boundaries_ccpp(low_boundaries, high_boundaries, units, errmsg, errflg)

   ! provide spectral boundaries of each longwave band

   real(kind_phys), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
   character(*), intent(in) :: units ! requested units
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   character(len=*), parameter :: sub = 'get_lw_spectral_boundaries_ccpp'
   !----------------------------------------------------------------------------

   ! Set error variables
   errmsg = ''
   errflg = 0

   if (.not. wavenumber_boundaries_set) then
      write(errmsg,'(a,a)') sub, ': ERROR, wavenumber boundaries not set.'
   end if

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries  = wavenumber_low_longwave
      high_boundaries = wavenumber_high_longwave
   case('m','meter','meters')
      low_boundaries  = 1.e-2_kind_phys/wavenumber_high_longwave
      high_boundaries = 1.e-2_kind_phys/wavenumber_low_longwave
   case('nm','nanometer','nanometers')
      low_boundaries  = 1.e7_kind_phys/wavenumber_high_longwave
      high_boundaries = 1.e7_kind_phys/wavenumber_low_longwave
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries  = 1.e4_kind_phys/wavenumber_high_longwave
      high_boundaries = 1.e4_kind_phys/wavenumber_low_longwave
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._kind_phys/wavenumber_high_longwave
      high_boundaries = 1._kind_phys/wavenumber_low_longwave
   case default
      write(errmsg, '(a,a,a)') sub, ': ERROR, requested spectral units not recognized: ', units
      errflg = 1
   end select

end subroutine get_lw_spectral_boundaries_ccpp

!=========================================================================================

subroutine get_mu_lambda_weights_ccpp(nmu, nlambda, g_mu, g_lambda, lamc, pgam, &
                mu_wgts, lambda_wgts, errmsg, errflg)
  integer,         intent(in) :: nmu
  integer,         intent(in) :: nlambda
  real(kind_phys), intent(in) :: g_mu(:)
  real(kind_phys), intent(in) :: g_lambda(:,:)
  real(kind_phys), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(kind_phys), intent(in) :: pgam   ! prognosed value of mu for cloud
  ! Output interpolation weights. Caller is responsible for freeing these.
  type(interp_type), intent(out) :: mu_wgts
  type(interp_type), intent(out) :: lambda_wgts
  character(len=*),  intent(out) :: errmsg
  integer,           intent(out) :: errflg

  integer :: ilambda
  real(kind_phys) :: g_lambda_interp(nlambda)

  ! Set error variables
  errmsg = ''
  errflg = 0

  ! Make interpolation weights for mu.
  ! (Put pgam in a temporary array for this purpose.)
  call lininterp_init(g_mu, nmu, [pgam], 1, extrap_method_bndry, mu_wgts)

  ! Use mu weights to interpolate to a row in the lambda table.
  do ilambda = 1, nlambda
     call lininterp(g_lambda(:,ilambda), nmu, &
          g_lambda_interp(ilambda:ilambda), 1, mu_wgts)
  end do

  ! Make interpolation weights for lambda.
  call lininterp_init(g_lambda_interp, nlambda, [lamc], 1, &
       extrap_method_bndry, lambda_wgts)

end subroutine get_mu_lambda_weights_ccpp

!=========================================================================================

end module radiation_utils
