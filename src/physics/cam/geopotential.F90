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

use ppgrid,                    only: pcols
use constituents,              only: pcnst, cnst_get_ind
use ccpp_constituent_prop_mod, only: ccpp_const_props      !CCPP constituent properties array (CAM version)
use geopotential_temp,         only: geopotential_temp_run !CCPP version
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
    logical  :: lagrang                 ! Lagrangian vertical coordinate flag
    integer  :: ixq                     ! state constituent array index for water vapor
    integer  :: i,k,idx                 ! Lon, level indices, water species index
    real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(ncol)               ! off-diagonal element
    real(r8) :: rog(ncol,pver)          ! Rair / gravit
    real(r8) :: tv                      ! virtual temperature
    real(r8) :: tvfac                   ! Tv/T
    real(r8) :: qfac(ncol,pver)         ! factor to convert from wet to dry mixing ratio
    real(r8) :: sum_dry_mixing_ratio(ncol,pver)! sum of dry water mixing ratios

    !CCPP-required variables (not used):
    integer            :: errflg
    character(len=512) :: errmsg

!
!-----------------------------------------------------------------------
!
    !Determine index for water vapor mass mixing ratio
    call cnst_get_ind('Q', ixq)

    !
    ! original code for backwards compatability with FV and EUL
    !
    if (.not.(dycore_is('MPAS') .or. dycore_is('SE'))) then

      !dry air gas constant over gravity
      rog(:ncol,:) = rair(:ncol,:) / gravit

      ! The surface height is zero by definition.
      do i = 1,ncol
        zi(i,pverp) = 0.0_r8
      end do

      ! Compute zi, zm from bottom up.
      ! Note, zi(i,k) is the interface above zm(i,k)
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
          tvfac   = 1._r8 + zvir(i,k) * q(i,k,ixq)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
        end do
      end do
    else !Using MPAS or SE dycore

      !Determine vertical coordinate type,
      !NOTE: Currently the FV (LR) or FV3 dycores
      !      do not allow for condensate loading,
      !      so for now 'lagrang' will always be FALSE.
      if ((dycore_is('LR') .or. dycore_is('FV3'))) then
        lagrang = .true.
      else
        lagrang = .false.
      end if

      !Use CCPP version of geopotential_t:
      call geopotential_temp_run(pver, lagrang, pver, 1, pverp, 1,  &
        pcnst, piln(1:ncol,:), pint(1:ncol,:), pmid(1:ncol,:),      &
        pdel(1:ncol,:), rpdel(1:ncol,:), t(1:ncol,:),               &
        q(1:ncol,:,ixq), q(1:ncol,:,:), ccpp_const_props,           &
        rair(1:ncol,:), gravit, zvir(1:ncol,:), zi(1:ncol,:),       &
        zm(1:ncol,:), ncol, errflg, errmsg)

    end if
  end subroutine geopotential_t
end module geopotential
