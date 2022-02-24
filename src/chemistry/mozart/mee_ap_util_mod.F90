!------------------------------------------------------------------------------
! Utility routines for medium energy electron ionization base on Ap
! geomagnetic activity index or observed electron fluxes
!------------------------------------------------------------------------------
module mee_ap_util_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: pi => shr_const_pi
  use mee_fluxes, only: mee_fluxes_nenergy, mee_fluxes_energy, mee_fluxes_denergy
  use mee_fluxes, only: mee_fluxes_active, mee_fluxes_extract
  use cam_abortutils, only : endrun

  implicit none

  private
  public :: mee_ap_iprs
  public :: mee_ap_init
  public :: mee_ap_final
  public :: mee_ap_error
  public :: mee_ap_noerror

  integer, parameter :: mee_ap_error=-1
  integer, parameter :: mee_ap_noerror=0

  integer :: nbins=0

  real(r8),pointer :: energies(:) => null()  ! energy bins
  real(r8),pointer :: denergies(:) => null() ! width of each energy bin
  real(r8) :: solid_angle_factor = -huge(1._r8)
  real(r8),allocatable :: fang_coefs(:,:)

contains

  !-----------------------------------------------------------------------------
  ! Sets up energy spectrum grid for medium energy electrons in earth's
  ! radiation belt
  !-----------------------------------------------------------------------------
  subroutine mee_ap_init(loss_cone_angle, status)

    real(r8), intent(in) :: loss_cone_angle ! Bounce Loss Cone angle (degrees)
    integer, intent(out) :: status          ! error status

    integer :: ierr
    character(len=*), parameter :: subname = 'mee_ap_init: '

    status = mee_ap_noerror
    if ( loss_cone_angle<0._r8 .or. loss_cone_angle>90._r8 ) then
       status = mee_ap_error
       return
    endif

    ! The area of the BLC in sr is 2pi(1-cos(BLC))
    solid_angle_factor = 2._r8*pi*(1._r8-cos(loss_cone_angle*pi/180._r8))

    if (mee_fluxes_active) then
       nbins=mee_fluxes_nenergy
       energies=>mee_fluxes_energy
       denergies=>mee_fluxes_denergy
    else
       nbins=100
       allocate(energies(nbins), stat=ierr)
       if (ierr/=0) call endrun(subname//'not able to allocate energies')
       allocate(denergies(nbins), stat=ierr)
       if (ierr/=0) call endrun(subname//'not able to allocate denergies')
       call gen_energy_grid(nbins, energies, denergies)
    endif

    call init_fang_coefs()

  end subroutine mee_ap_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ap_final()

    deallocate(fang_coefs)

    if (.not.mee_fluxes_active) then
       deallocate(energies)
       deallocate(denergies)
    endif
    nullify(energies)
    nullify(denergies)

 end subroutine mee_ap_final

  !-----------------------------------------------------------------------------
  ! Computes ion pair production rates base on Ap geomagnetic activity index
  ! or observed electron fluxes
  !-----------------------------------------------------------------------------
  subroutine mee_ap_iprs( ncols, nlyrs, airden, scaleh, Ap, ionpairs, status, maglat, lshell)

    integer ,intent(in) :: ncols, nlyrs
    real(r8),intent(in) :: airden(ncols,nlyrs)     ! g/cm3
    real(r8),intent(in) :: scaleh(ncols,nlyrs)     ! cm
    real(r8),intent(in) :: Ap                      ! geomagnetic activity index
    real(r8),intent(out) :: ionpairs(ncols,nlyrs)  ! pairs/cm3/sec
    integer, intent(out) :: status                 ! error status
    real(r8),intent(in), optional :: maglat(ncols) ! magnetic latitude (radians)
    real(r8),intent(in), optional :: lshell(ncols) ! magnetosphere L-Shells

    integer :: i,k
    real(r8) :: flux_sd(nbins)
    logical  :: valid(nbins)
    real(r8) :: flux(nbins)
    real(r8) :: ipr(nbins,nlyrs)
    real(r8) :: l_shells(ncols)

    status = 0

    if (present(lshell)) then
       l_shells(:) = lshell(:)
    elseif (present(maglat)) then
       ! get L-shell values corresponeding to the column geo-mag latitudes
       l_shells = maglat2lshell( maglat )
    else
       ionpairs(:,:) = -huge(1._r8)
       status = mee_ap_error
       return
    endif

    ionpairs(:,:) = 0._r8

    do i = 1,ncols

       if (  l_shells(i) >= 2._r8 .and. l_shells(i) <= 10._r8 ) then

          valid(:) = .false.
          flux_sd(:) = 0._r8

          if (mee_fluxes_active) then
             ! use prescribed top-of-atmosphere fluxes
             call mee_fluxes_extract( l_shells(i), flux_sd, valid )
          end if

          where ( (.not.valid(:)) .and. (energies(:)>=30._r8) .and. (energies(:)<=1000._r8) )
             ! calculate the top of the atmosphere energetic electron energy spectrum
             ! van de Kamp is per steradian (electrons / (cm2 sr s keV))
             flux_sd(:) = FluxSpectrum(energies(:), l_shells(i), Ap)
          end where

          ! assume flux is isotropic inside a nominal bounce loss cone (BLC) angle.
          ! The area of the BLC in sr is 2pi(1-cosd(BLC))
          flux(:) = solid_angle_factor*flux_sd(:)

          ! calculate the IPR as a function of height and energy
          ipr(:,:) = iprmono(energies, flux, airden(i,:), scaleh(i,:))

          ! integrate across the energy range to get total IPR
          do k=1,nlyrs
             ionpairs(i,k) = sum(ipr(:,k)*denergies(:))
          end do

       end if

    end do

  end subroutine mee_ap_iprs

  !------------------------------------------------------------------------------
  ! Electron fluxes for a specific L-shell and Ap
  !
  ! Based on:
  !
  ! van de Kamp, M., Seppala, A., Clilverd, M. A., Rodger, C. J., Verronen, P. T., and Whittaker, I. C. (2016),
  ! A model providing long‐term datasets of energetic electron precipitation during geomagnetic storms,
  ! J. Geophys. Res. Atmos., 121, 12,520– 12,540,
  ! [doi:10.1002/2015JD024212](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD024212)
  !------------------------------------------------------------------------------
  pure function FluxSpectrum( en, lshell, Ap ) result(flux)

    real(r8), intent(in) :: en(:)  ! electron energy bins (keV)
    real(r8), intent(in) :: lshell ! magnetosphere L-Shell number
    real(r8), intent(in) :: Ap     ! geomagnetic activity index

    real(r8) :: flux(size(en))

    real(r8) :: lpp
    real(r8) :: Spp
    real(r8) :: A
    real(r8) :: b
    real(r8) :: c
    real(r8) :: s
    real(r8) :: d
    real(r8) :: F30
    real(r8) :: E
    real(r8) :: bk
    real(r8) :: sk
    real(r8) :: k
    real(r8) :: x

    lpp = -0.7430_r8*log(Ap) + 6.5257_r8
    Spp = lshell - lpp

    ! vdK2016 eqn.(8)

    A = 8.2091_r8*Ap**(0.16255_r8)
    b = 1.3754_r8*Ap**(0.33042_r8)
    c = 0.13334_r8*Ap**(0.42616_r8)
    s = 2.2833_r8*Ap**(-0.22990_r8)
    d = 2.7563e-4_r8*Ap**(2.6116_r8)

    F30 = exp(A) / (exp(-b*(Spp-s)) + exp(c*(Spp-s)) + d)

    ! vdK2016 eqn.(9)

    E = 3.3777_r8*Ap**(-1.7038_r8) + 0.15_r8
    bk = 3.7632_r8*Ap**(-0.16034_r8)
    sk = 12.184_r8*Ap**(-0.30111_r8)

    k = -1.0_r8 / (E*exp(-bk*Spp) + 0.30450_r8*cosh(0.20098_r8*(Spp-sk))) - 1._r8

    x=k+1
    c = F30*x / ((1e3_r8**x) - (30._r8**x))
    flux(:) = c * (en(:)**k)

  end function FluxSpectrum

  !------------------------------------------------------------------------------
  ! Calculate how energy from top of atmosphere is deposited in rest of atmosphere
  !
  ! The function takes the energy spectrum at the top of the atmosphere and
  ! calculates how that energy is deposited in the atmosphere using the parameterization
  ! described in [Fang et al., (2010)](https://opensky.ucar.edu/islandora/object/articles:10653)
  !
  ! Fang, X., C. E. Randall, D. Lummerzheim, W. Wang, G. Lu, S. C. Solomon,
  ! and R. A. Frahm (2010), Parameterization of monoenergetic electron impact
  ! ionization, Geophys. Res. Lett., 37, L22106, [doi:10.1029/2010GL045406.]
  ! (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2010GL045406)
  !
  ! Application of the new parameterization requires the following steps:
  !
  ! 1. Fang coefficients using equation (5) and Table 1. are calculated at initialization time
  ! 2. Calculate the y values throughout the atmosphere using equation (1).
  ! 3. Calculate the normalized energy dissipation f values using equation (4).
  ! 4. Obtain the altitude profile of qtot by substituting the f values into equation (3).
  !
  !------------------------------------------------------------------------------
  function iprmono(e, flux, rho, scaleh) result(ipr)
    real(r8), intent(in) :: e(:)      ! electron energy bins (keV)
    real(r8), intent(in) :: flux(:)   ! top of atmos electron fluxes (electrons / (cm2 sr s keV))
    real(r8), intent(in) :: rho(:)    ! density of atmos (g/cm3)
    real(r8), intent(in) :: scaleh(:) ! scale  height (cm)

    real(r8) :: ipr(size(e),size(rho))
    integer :: nenergies, n
    integer :: nlayers, k

    real(r8) :: y(size(rho))
    real(r8) :: f(size(rho))
    real(r8) :: Qmono

    ! assign constants
    real(r8), parameter :: epsilon = 0.035_r8 ! keV energy loss per ion pair produced

    ipr = 0._r8
    nenergies = size(e)
    nlayers = size(rho)

    do n = 1,nenergies

       ! step 1. (eq. 1)
       y(:) = (2/e(n))*(rho(:)*scaleh(:)/6.0e-6_r8)**0.7_r8

       do k = 1,nlayers
          f(k) = fang(y(k), n)
       end do

       ! calculate ipr (qtot) using eq. (3) for a specified flux at ea. energy
       Qmono = flux(n)*e(n) ! (keV cm−2 s−1)
       ipr(n,:) = f(:)*Qmono/(epsilon*scaleh(:))
    end do

  contains

    pure function fang(y,n) result(f)
      real(r8), intent(in) :: y
      integer,  intent(in) :: n ! energy ndx

      real(r8) :: f
      ! Input:
      ! y - normalized atmospheric column mass as a function of vertical location (z)
      ! Emono - is incident electron energy (keV)
      ! Output:
      ! f - quanity calculated by eqn. (4)


      ! eq. (4) - Normalized energy deposition
      f = fang_coefs(1,n) * (y**fang_coefs(2,n)) * exp(-fang_coefs(3,n)*y**fang_coefs(4,n)) &
        + fang_coefs(5,n) * (y**fang_coefs(6,n)) * exp(-fang_coefs(7,n)*y**fang_coefs(8,n))

    end function fang

  end function iprmono

  !------------------------------------------------------------------------------
  ! initializes the coeffs used in the fang function in iprmono
  !------------------------------------------------------------------------------
  subroutine init_fang_coefs

    integer :: n,i, ierr

    real(r8) :: lne, lne2, lne3
    ! Table 1. of Fang et al., 2010
    ! (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2010GL045406)
    real(r8), parameter :: p1(8,4) = reshape( &
         (/ 1.24616E+0_r8,  1.45903E+0_r8, -2.42269E-1_r8,  5.95459E-2_r8, &
            2.23976E+0_r8, -4.22918E-7_r8,  1.36458E-2_r8,  2.53332E-3_r8, &
            1.41754E+0_r8,  1.44597E-1_r8,  1.70433E-2_r8,  6.39717E-4_r8, &
            2.48775E-1_r8, -1.50890E-1_r8,  6.30894E-9_r8,  1.23707E-3_r8, &
           -4.65119E-1_r8, -1.05081E-1_r8, -8.95701E-2_r8,  1.22450E-2_r8, &
            3.86019E-1_r8,  1.75430E-3_r8, -7.42960E-4_r8,  4.60881E-4_r8, &
           -6.45454E-1_r8,  8.49555E-4_r8, -4.28581E-2_r8, -2.99302E-3_r8, &
            9.48930E-1_r8,  1.97385E-1_r8, -2.50660E-3_r8, -2.06938E-3_r8 /), &
         shape=(/8,4/),order=(/2,1/))

    allocate(fang_coefs(8,nbins), stat=ierr)
    if (ierr/=0) call endrun('init_fang_coefs: not able to allocate fang_coefs')

    do n = 1,nbins
       ! terms in eq. (5)
       lne = log(energies(n))
       lne2 = lne*lne
       lne3 = lne*lne2

       ! step 2. calculate the C array in (5)
       do i = 1,8
          fang_coefs(i,n) = exp(p1(i,1) + p1(i,2)*lne + p1(i,3)*lne2 + p1(i,4)*lne3)
       end do

    end do
  end subroutine init_fang_coefs

  !------------------------------------------------------------------------------
  ! Generate a grid of energy bins for the flux spectrum.
  ! The energy range of the spectrum is 30-1000 keV,
  ! with nbins of logarithmically spaced grid points.
  !------------------------------------------------------------------------------
  subroutine gen_energy_grid(nbins, energies, deltas)
    integer, intent(in) :: nbins
    real(r8),intent(out) :: energies(nbins)
    real(r8),intent(out) :: deltas(nbins)

    integer :: i
    real(r8) :: low,med,hig

    ! for energy bins ranging from 30 keV to 1000 keV
    real(r8), parameter :: e1 = log10(30._r8)
    real(r8), parameter :: e2 = log10(1000._r8)

    do i = 1,nbins
       low = e1 + (e2-e1)*(i-1.0_r8)/nbins
       med = e1 + (e2-e1)*(i-0.5_r8)/nbins
       hig = e1 + (e2-e1)*(i)/nbins

       energies(i) = 10**med
       deltas(i) = (10**hig)-(10**low)
    end do

  end subroutine gen_energy_grid

  !------------------------------------------------------------------------------
  ! returns L-Shell number for a given magnetic latitude (radians)
  !------------------------------------------------------------------------------
  pure elemental function maglat2lshell( rmaglat ) result(lshell)
    real(r8), intent(in) :: rmaglat ! mag latitude in radians
    real(r8) :: lshell  ! magnetosphere L-Shell number

    ! lshell = r/(cos(rmaglat)**2) (https://en.wikipedia.org/wiki/L-shell)
    ! where r is the radial distance (in planetary radii) to a point on the line.
    ! r = 1.01 corresponds to an altitude in the lower mesosphere (~64 km)
    ! where medium-energy electrons typically deposit their energy.
    lshell = 1.01_r8/(cos(rmaglat)**2)

  end function maglat2lshell

end module mee_ap_util_mod
