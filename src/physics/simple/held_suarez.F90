module held_suarez
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  ! 
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  public :: held_suarez_1994_init
  public :: held_suarez_1994

  !!
  !! Forcing parameters
  !!
  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

  !!
  !! Model constants, reset in init call
  !!
  real(r8)              :: cappa = 2.0_r8 / 7.0_r8  ! R/Cp
  real(r8)              :: cpair = 1004.0_r8        ! specific heat of dry air (J/K/kg)
  real(r8)              :: psurf_ref = 0.0_r8       ! Surface pressure
  ! pref_mid_norm are layer midpoints normalized by surface pressure ('eta' coordinate)
  real(r8), allocatable :: pref_mid_norm(:)
  integer               :: pver                     ! Num vertical levels



!======================================================================= 
contains
!======================================================================= 

  subroutine held_suarez_1994_init(cappa_in, cpair_in, psurf_ref_in, pref_mid_norm_in)
    !! Dummy arguments
    real(r8), intent(in) :: cappa_in
    real(r8), intent(in) :: cpair_in
    real(r8), intent(in) :: psurf_ref_in
    real(r8), intent(in) :: pref_mid_norm_in(:)

    pver = size(pref_mid_norm_in)
    allocate(pref_mid_norm(pver))
    cappa         = cappa_in
    cpair         = cpair_in
    psurf_ref     = psurf_ref_in
    pref_mid_norm = pref_mid_norm_in

  end subroutine held_suarez_1994_init

  subroutine held_suarez_1994(pcols, ncol, clat, pmid, &
       u, v, t, du, dv, s)

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pcols            ! Size of column dimension
    integer,  intent(in)  :: ncol             ! Num active columns
    real(r8), intent(in)  :: clat(pcols)      ! latitudes(radians) for columns
    real(r8), intent(in)  :: pmid(pcols,pver) ! mid-point pressure
    real(r8), intent(in)  :: u(pcols,pver)    ! Zonal wind (m/s)
    real(r8), intent(in)  :: v(pcols,pver)    ! Meridional wind (m/s)
    real(r8), intent(in)  :: t(pcols,pver)    ! Temperature (K)
                                              !
                                              ! Output arguments
                                              !
    real(r8), intent(out) :: du(pcols,pver)   ! Zonal wind tend
    real(r8), intent(out) :: dv(pcols,pver)   ! Meridional wind tend
    real(r8), intent(out) :: s(pcols,pver)    ! Heating rate
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(r8) :: kv            ! 1./efolding_time (normalized) for wind
    real(r8) :: kt            ! 1./efolding_time for temperature diss.
    real(r8) :: trefa         ! "radiative equilibrium" T
    real(r8) :: trefc         ! used in calc of "radiative equilibrium" T
    real(r8) :: cossq(ncol)   ! coslat**2
    real(r8) :: cossqsq(ncol) ! coslat**4
    real(r8) :: sinsq(ncol)   ! sinlat**2
    real(r8) :: coslat(ncol)  ! cosine(latitude)
    !
    !-----------------------------------------------------------------------
    !

    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
    end do

    !
    !-----------------------------------------------------------------------
    !
    ! Held/Suarez IDEALIZED physics algorithm:
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the
    !   intercomparison of the dynamical cores of atmospheric general
    !   circulation models.
    !   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
    !
    !-----------------------------------------------------------------------
    !
    ! Compute idealized radiative heating rates (as dry static energy)
    !
    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        do i = 1, ncol
          kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
          trefc   = 315._r8 - (60._r8 * sinsq(i))
          trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*kt*cpair
        end do
      else
        do i = 1, ncol
          trefc   = 315._r8 - 60._r8*sinsq(i)
          trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*ka*cpair
        end do
      end if
    end do
    !
    ! Add diffusion near the surface for the wind fields
    !
    do k = 1, pver
      do i = 1, pcols
        du(i,k) = 0._r8
        dv(i,k) = 0._r8
      end do
    end do

    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
        do i = 1, ncol
          du(i,k) = -kv*u(i,k)
          dv(i,k) = -kv*v(i,k)
        end do
      end if
    end do

  end subroutine held_suarez_1994

end module held_suarez
