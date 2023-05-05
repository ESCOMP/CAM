
module geopotential

!---------------------------------------------------------------------------------
! Compute geopotential from temperature or
! compute geopotential and temperature from dry static energy.
!
! The hydrostatic matrix elements must be consistent with the dynamics algorithm.
! The diagonal element is the itegration weight from interface k+1 to midpoint k.
! The offdiagonal element is the weight between interfaces.
! 
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver, pverp
  use dycore,       only: dycore_is

  implicit none
  private
  save

  public geopotential_t

contains

!===============================================================================
  subroutine geopotential_t(                                 &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------

use ppgrid, only : pcols
use air_composition, only: thermodynamic_active_species_num,thermodynamic_active_species_idx
!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: ncol                  ! Number of longitudes

    real(r8), intent(in) :: piln (:,:)    ! (pcols,pverp) - Log interface pressures
    real(r8), intent(in) :: pmln (:,:)    ! (pcols,pver)  - Log midpoint pressures
    real(r8), intent(in) :: pint (:,:)    ! (pcols,pverp) - Interface pressures
    real(r8), intent(in) :: pmid (:,:)    ! (pcols,pver)  - Midpoint pressures
    real(r8), intent(in) :: pdel (:,:)    ! (pcols,pver)  - layer thickness
    real(r8), intent(in) :: rpdel(:,:)    ! (pcols,pver)  - inverse of layer thickness
    real(r8), intent(in) :: t    (:,:)    ! (pcols,pver)  - temperature
    real(r8), intent(in) :: q    (:,:,:)  ! (pcols,pver,:)- tracers (moist mixing ratios)
    real(r8), intent(in) :: rair (:,:)    ! (pcols,pver)  - Gas constant for dry air
    real(r8), intent(in) :: gravit        !               - Acceleration of gravity
    real(r8), intent(in) :: zvir (:,:)    ! (pcols,pver)  - rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: zi(:,:)      ! (pcols,pverp) - Height above surface at interfaces
    real(r8), intent(out) :: zm(:,:)      ! (pcols,pver)  - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    integer  :: i,k,idx                 ! Lon, level indices, water species index
    real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(ncol)               ! off-diagonal element
    real(r8) :: rog(ncol,pver)          ! Rair / gravit
    real(r8) :: tv                      ! virtual temperature
    real(r8) :: tvfac                   ! Tv/T
    real(r8) :: qfac(ncol,pver)         ! factor to convert from wet to dry mixing ratio
    real(r8) :: sum_dry_mixing_ratio(ncol,pver)! sum of dry water mixing ratios

!
!-----------------------------------------------------------------------
!
    rog(:ncol,:) = rair(:ncol,:) / gravit

! The surface height is zero by definition.

    do i = 1,ncol
       zi(i,pverp) = 0.0_r8
    end do

    ! Compute zi, zm from bottom up. 
    ! Note, zi(i,k) is the interface above zm(i,k)

    !
    ! original code for backwards compatability with FV and EUL
    !
    if (.not.(dycore_is('MPAS') .or. dycore_is('SE'))) then
      do k = pver, 1, -1
        
        ! First set hydrostatic elements consistent with dynamics
        
        if ((dycore_is('LR') .or. dycore_is('FV3'))) then
          do i = 1,ncol
            hkl(i) = piln(i,k+1) - piln(i,k)
            hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          end do
        else
          do i = 1,ncol
            hkl(i) = pdel(i,k) / pmid(i,k)
            hkk(i) = 0.5_r8 * hkl(i)
          end do
        end if
        
        ! Now compute tv, zm, zi
        
        do i = 1,ncol
          tvfac   = 1._r8 + zvir(i,k) * q(i,k,1)
          tv      = t(i,k) * tvfac
          
          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
        end do
      end do
    else
      !
      ! For the computation of generalized virtual temperature (equation 16
      ! in Lauritzen et al. (2018);  https://doi.org/10.1029/2017MS001257)
      !
      ! Compute factor for converting wet to dry mixing ratio (eq.7)
      !
      qfac = 1.0_r8
      do idx = 1,thermodynamic_active_species_num
        do k = 1,pver
          do i = 1,ncol
            qfac(i,k) = qfac(i,k)-q(i,k,thermodynamic_active_species_idx(idx))
          end do
        end do
      end do
      qfac = 1.0_r8/qfac
      
      ! Compute sum of dry water mixing ratios
      sum_dry_mixing_ratio = 1.0_r8
      do idx = 1,thermodynamic_active_species_num
        do k = 1,pver
          do i = 1,ncol
            sum_dry_mixing_ratio(i,k) = sum_dry_mixing_ratio(i,k)&
                 +q(i,k,thermodynamic_active_species_idx(idx))*qfac(i,k)
          end do
        end do
      end do
      sum_dry_mixing_ratio(:,:) = 1.0_r8/sum_dry_mixing_ratio(:,:)
      ! Compute zi, zm from bottom up. 
      ! Note, zi(i,k) is the interface above zm(i,k)
      do k = pver, 1, -1
        
        ! First set hydrostatic elements consistent with dynamics
 
        !
        ! the outcommented code is left for when/if we will support
        ! FV3 and/or FV with condensate loading
        !
       
!        if ((dycore_is('LR') .or. dycore_is('FV3'))) then
!          do i = 1,ncol
!            hkl(i) = piln(i,k+1) - piln(i,k)
!            hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
!          end do
!        else!MPAS, SE or EUL
          !
          ! For SE   : pmid = 0.5*(pint(k+1)+pint(k))
          ! For MPAS : pmid is computed from theta_m, rhodry, etc.
          !
          do i = 1,ncol
            hkl(i) = pdel(i,k) / pmid(i,k)
            hkk(i) = 0.5_r8 * hkl(i)
          end do
!        end if
        
        ! Now compute tv, zm, zi
        
        do i = 1,ncol
          tvfac   = (1._r8 + (zvir(i,k)+1.0_r8) * q(i,k,1)*qfac(i,k))*sum_dry_mixing_ratio(i,k)
          tv      = t(i,k) * tvfac
          
          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
        end do
      end do
    end if
    return
  end subroutine geopotential_t
end module geopotential
