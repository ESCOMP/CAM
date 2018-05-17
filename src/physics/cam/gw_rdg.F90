module gw_rdg

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
use shr_const_mod, only: pii => shr_const_pi
use gw_utils, only: r8
use gw_common, only: pver
use coords_1d, only: Coords1D
use spmd_utils,only: masterproc
use cam_abortutils, only: endrun


implicit none
private
save

! Public interface
public :: gw_rdg_readnl
public :: gw_rdg_src
public :: gw_rdg_belowpeak
public :: gw_rdg_break_trap
public :: gw_rdg_do_vdiff

! Tunable Parameters
!--------------------
logical            :: do_divstream

!===========================================
! Parameters for DS2017 (do_divstream=.T.)
!===========================================
! Amplification factor - 1.0 for
! high-drag/windstorm regime
real(r8), protected :: C_BetaMax_DS

! Max Ratio Fr2:Fr1 - 1.0
real(r8), protected :: C_GammaMax

! Normalized limits  for Fr2(Frx) function
real(r8), protected :: Frx0
real(r8), protected :: Frx1


!===========================================
! Parameters for SM2000
!===========================================
! Amplification factor - 1.0 for
! high-drag/windstorm regime
real(r8), protected :: C_BetaMax_SM



! NOTE: Critical inverse Froude number Fr_c is 
! 1./(SQRT(2.)~0.707 in SM2000
! (should be <= 1)
real(r8), protected :: Fr_c

logical, protected :: gw_rdg_do_vdiff=.true.

logical :: do_smooth_regimes
logical :: do_adjust_tauoro
logical :: do_backward_compat


! Limiters (min/max values)
! min surface displacement height for orographic waves (m)
real(r8), protected :: orohmin
! min wind speed for orographic waves
real(r8), protected :: orovmin
! min stratification allowing wave behavior
real(r8), protected :: orostratmin
! min stratification allowing wave behavior
real(r8), protected :: orom2min

!==========================================================================
contains
!==========================================================================

subroutine gw_rdg_readnl(nlfile)
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, mpi_logical

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gw_rdg_readnl'

  logical ::  gw_rdg_do_divstream, gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
              gw_rdg_do_backward_compat

  
  real(r8) :: gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, &
              gw_rdg_Frx0, gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, &
              gw_rdg_orohmin, gw_rdg_orovmin, gw_rdg_orostratmin, gw_rdg_orom2min 

  namelist /gw_rdg_nl/ gw_rdg_do_divstream, gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, &
                       gw_rdg_Frx0, gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, &
                       gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
                       gw_rdg_do_backward_compat, gw_rdg_orohmin, gw_rdg_orovmin, &
                       gw_rdg_orostratmin, gw_rdg_orom2min, gw_rdg_do_vdiff

  !----------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'gw_rdg_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, gw_rdg_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(sub // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

     ! Set the local variables
     do_divstream        = gw_rdg_do_divstream 
     C_BetaMax_DS        = gw_rdg_C_BetaMax_DS
     C_GammaMax          = gw_rdg_C_GammaMax
     Frx0                = gw_rdg_Frx0
     Frx1                = gw_rdg_Frx1
     C_BetaMax_SM        = gw_rdg_C_BetaMax_SM
     Fr_c                = gw_rdg_Fr_c
     do_smooth_regimes   = gw_rdg_do_smooth_regimes
     do_adjust_tauoro    = gw_rdg_do_adjust_tauoro
     do_backward_compat  = gw_rdg_do_backward_compat
     orohmin             = gw_rdg_orohmin
     orovmin             = gw_rdg_orovmin
     orostratmin         = gw_rdg_orostratmin
     orom2min            = gw_rdg_orom2min
  end if

  ! Broadcast the local variables

  call mpi_bcast(do_divstream, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_divstream")
  call mpi_bcast(do_smooth_regimes, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_smooth_regimes")
  call mpi_bcast(do_adjust_tauoro, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_adjust_tauoro")
  call mpi_bcast(do_backward_compat, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_backward_compat")

  call mpi_bcast(C_BetaMax_DS, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: C_BetaMax_DS")
  call mpi_bcast(C_GammaMax, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: C_GammaMax")
  call mpi_bcast(Frx0, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: Frx0")
  call mpi_bcast(Frx1, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: Frx1")
  call mpi_bcast(C_BetaMax_SM, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: C_BetaMax_SM")
  call mpi_bcast(Fr_c, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: Fr_c")
  call mpi_bcast(orohmin, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: orohmin")
  call mpi_bcast(orovmin, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: orovmin")
  call mpi_bcast(orostratmin, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: orostratmin")
  call mpi_bcast(orom2min, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: orom2min")

  call mpi_bcast(gw_rdg_do_vdiff, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_vdiff")

  if (Fr_c > 1.0_r8) call endrun(sub//": FATAL: Fr_c must be <= 1")

end subroutine gw_rdg_readnl


subroutine gw_rdg_src(ncol, band, p, &
     u, v, t, mxdis, angxy, anixy, kwvrdg, iso, zi, nm, &
     src_level, tend_level, bwv_level ,tlb_level , tau, ubm, ubi, xv, yv,  & 
     ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, c)
  use gw_common, only: rair, GWBand
  use gw_utils, only:  dot_2d, midpoint_interp
  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol

  ! Band to emit orographic waves in.
  ! Regardless, we will only ever emit into l = 0.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p


  ! Midpoint zonal/meridional winds. ( m s-1)
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Midpoint temperatures. (K)
  real(r8), intent(in) :: t(ncol,pver)
  ! Height estimate for ridge (m) [anisotropic orography].
  real(r8), intent(in) :: mxdis(ncol)
  ! Angle of ridge axis w/resp to north (degrees) [anisotropic orography].
  real(r8), intent(in) :: angxy(ncol)
  ! Anisotropy parameter [anisotropic orography].
  real(r8), intent(in) :: anixy(ncol)
  ! horiz wavenumber [anisotropic orography].
  real(r8), intent(in) :: kwvrdg(ncol)
  ! Isotropic source flag [anisotropic orography].
  integer, intent(in)  :: iso(ncol)
  ! Interface altitudes above ground (m).
  real(r8), intent(in) :: zi(ncol,pver+1)
  ! Midpoint Brunt-Vaisalla frequencies (s-1).
  real(r8), intent(in) :: nm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)
  integer, intent(out) :: bwv_level(ncol),tlb_level(ncol)

  ! Averages over source region.
  real(r8), intent(out) :: nsrc(ncol) ! B-V frequency.
  real(r8), intent(out) :: rsrc(ncol) ! Density.
  real(r8), intent(out) :: usrc(ncol) ! Zonal wind.
  real(r8), intent(out) :: vsrc(ncol) ! Meridional wind.
  real(r8), intent(out) :: ubmsrc(ncol) ! On-ridge wind.
  ! Top of low-level flow layer.
  real(r8), intent(out) :: tlb(ncol)
  ! Bottom of linear wave region.
  real(r8), intent(out) :: bwv(ncol)
  ! normalized wavenumber
  real(r8), intent(out) :: m2src(ncol)


  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)
  ! Froude numbers for flow/drag regimes
  real(r8), intent(out) :: Fr1(ncol), Fr2(ncol), Frx(ncol)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Surface streamline displacement height (2*sgh).
  real(r8) :: hdsp(ncol)

  ! Difference in interface pressure across source region.
  real(r8) :: dpsrc(ncol)
  ! Thickness of downslope wind region.
  real(r8) :: ddw(ncol)
  ! Thickness of linear wave region.
  real(r8) :: dwv(ncol)
  ! Wind speed in source region.
  real(r8) :: wmsrc(ncol)

  real(r8) :: ragl(ncol) 
  
!--------------------------------------------------------------------------
! Check that ngwav is equal to zero, otherwise end the job
!--------------------------------------------------------------------------
  if (band%ngwv /= 0) call endrun(' gw_rdg_src :: ERROR - band%ngwv must be zero and it is not')

!--------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the appropiate
! values of wind, stability, etc. for determining the wave source are
! averages over the depth of the atmosphere penterated by the typical
! mountain.
! Reduces to the bottom midpoint values when mxdis=0, such as over ocean.
!--------------------------------------------------------------------------

  hdsp      = mxdis ! no longer multipied by 2
  src_level = pver+1
  bwv_level = -1
  tlb_level = -1

  tau(:,0,:) = 0.0_r8

  ! Find depth of "source layer" for mountain waves
  ! i.e., between ground and mountain top
  do k = pver, 1, -1
     do i = 1, ncol
        ! Need to have h >= z(k+1) here or code will bomb when h=0.
        if ( (hdsp(i) >= zi(i,k+1)) .and. (hdsp(i) < zi(i,k))   ) then
           src_level(i) = k  
        end if
     end do
  end do

  rsrc = 0._r8
  usrc = 0._r8 
  vsrc = 0._r8
  nsrc = 0._r8
  do i = 1, ncol
      do k = pver, src_level(i), -1
           rsrc(i) = rsrc(i) + p%mid(i,k) / (rair*t(i,k))* p%del(i,k)
           usrc(i) = usrc(i) + u(i,k) * p%del(i,k)
           vsrc(i) = vsrc(i) + v(i,k) * p%del(i,k)
           nsrc(i) = nsrc(i) + nm(i,k)* p%del(i,k)
     end do
  end do


  do i = 1, ncol
     dpsrc(i) = p%ifc(i,pver+1) - p%ifc(i,src_level(i))
  end do

  rsrc = rsrc / dpsrc
  usrc = usrc / dpsrc
  vsrc = vsrc / dpsrc
  nsrc = nsrc / dpsrc

  wmsrc = sqrt( usrc**2 + vsrc**2 )


  ! Get the unit vector components
  ! Want agl=0 with U>0 to give xv=1

  ragl = angxy * pii/180._r8

  ! protect from wierd "bad" angles 
  ! that may occur if hdsp is zero
  where( hdsp <= orohmin )
     ragl = 0._r8
  end where

  yv   =-sin( ragl )
  xv   = cos( ragl )


  ! Kluge in possible "isotropic" obstacle.
  where( ( iso == 1 ) .and. (wmsrc > orovmin) )
       xv = usrc/wmsrc    
       yv = vsrc/wmsrc
  end where


  ! Project the local wind at midpoints into the on-ridge direction
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do
  ubmsrc = dot_2d(usrc , vsrc , xv, yv)

  ! Ensure on-ridge wind is positive at source level
  do k = 1, pver
     ubm(:,k) = sign( ubmsrc*0._r8+1._r8 , ubmsrc ) *  ubm(:,k)
  end do

                  ! Sean says just use 1._r8 as 
                  ! first argument
  xv  = sign( ubmsrc*0._r8+1._r8 , ubmsrc ) *  xv
  yv  = sign( ubmsrc*0._r8+1._r8 , ubmsrc ) *  yv

  ! Now make ubmsrc positive and protect
  ! against zero
  ubmsrc = abs(ubmsrc)
  ubmsrc = max( 0.01_r8 , ubmsrc )
  

  ! The minimum stratification allowing GW behavior
  ! should really depend on horizontal scale since
  !
  !      m^2 ~ (N/U)^2 - k^2
  !
  ! Should also think about parameterizing
  ! trapped lee-waves.  

  
  ! This needs to be made constistent with later
  ! treatment of nonhydrostatic effects.
  m2src = ( (nsrc/(ubmsrc+0.01_r8))**2 - kwvrdg**2 ) /((nsrc/(ubmsrc+0.01_r8))**2)


  !-------------------------------------------------------------
  ! Calculate provisional limits (in Z [m]) for 3 regimes. This
  ! will modified later if wave breaking or trapping are
  ! diagnosed
  !
  !                                            ^ 
  !                                            | *** linear propagation ***
  !  (H) -------- mountain top -------------   | *** or wave breaking  ****     
  !                                            | *** regimes  *************
  ! (BWV)------ bottom of linear waves ----    |
  !                    :                       |
  !                 *******                    |
  !                    :                       |
  ! (TLB)--- top of flow diversion layer---    '
  !                   :
  !        **** flow diversion *****  
  !                    :
  !============================================

  !============================================
  ! For Dividing streamline para (DS2017)
  !--------------------------------------------
  ! High-drag downslope wind regime exists
  ! between bottom of linear waves and top of
  ! flow diversion. Linear waves can only 
  ! attain vertical displacment of f1*U/N. So,
  ! bottom of linear waves is given by
  !
  !        BWV = H - Fr1*U/N 
  !
  ! Downslope wind layer begins at BWV and 
  ! extends below it until some maximum high
  ! drag obstacle height Fr2*U/N is attained
  ! (where Fr2 >= f1).  Below downslope wind
  ! there is flow diversion, so top of 
  ! diversion layer (TLB) is equivalent to
  ! bottom of downslope wind layer and is;
  !
  !       TLB = H - Fr2*U/N
  !
  !-----------------------------------------

  ! Critical inverse Froude number
  !-----------------------------------------------
  Fr1(:) = Fr_c * 1.00_r8
  Frx(:) = hdsp(:)*nsrc(:)/abs( ubmsrc(:) ) / Fr_c

  if ( do_divstream ) then
     !------------------------------------------------
     ! Calculate Fr2(Frx) for DS2017   
     !------------------------------------------------
     where(Frx <= Frx0)
          Fr2(:) = Fr1(:) + Fr1(:)* C_GammaMax * anixy(:)
     elsewhere((Frx > Frx0).and.(Frx <= Frx1) )
          Fr2(:) = Fr1(:) + Fr1(:)* C_GammaMax * anixy(:) &
                   * (Frx1 - Frx(:))/(Frx1-Frx0)    
     elsewhere(Frx > Frx1) 
          Fr2(:)=Fr1(:)
     endwhere
  else
  !------------------------------------------   
  ! Regime distinctions entirely carried by
  ! amplification of taudsw (next subr)
  !------------------------------------------
     Fr2(:)=Fr1(:)
  end if   


  
  where( m2src > orom2min ) 
     ddw  = Fr2 * ( abs(ubmsrc) )/nsrc
  elsewhere
     ddw  = 0._r8
  endwhere


  ! If TLB is less than zero then obstacle is not
  ! high enough to produce an low-level diversion layer
  tlb = mxdis - ddw
  where( tlb < 0._r8)
     tlb = 0._r8
  endwhere
  do k = pver, pver/2, -1
     do i = 1, ncol
         if ( (tlb(i) > zi(i,k+1)) .and. (tlb(i) <= zi(i,k))   ) then
           tlb_level(i) = k
        end if
     end do
  end do


  ! Find *BOTTOM* of linear wave layer (BWV)
  !where ( nsrc > orostratmin )
  where( m2src > orom2min ) 
      dwv  = Fr1 * ( abs(ubmsrc) )/nsrc
  elsewhere
     dwv  = -9.999e9_r8 ! if weak strat - no waves
  endwhere

  bwv = mxdis - dwv
  where(( bwv < 0._r8) .or. (dwv < 0._r8) )
     bwv = 0._r8
  endwhere
  do k = pver,1, -1
     do i = 1, ncol
        if ( (bwv(i) > zi(i,k+1)) .and. (bwv(i) <= zi(i,k))   ) then
           bwv_level(i) = k+1
        end if
     end do
  end do



  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)
  ubi(:,2:pver) = midpoint_interp(ubm)
  ubi(:,pver+1) = ubm(:,pver)

  ! Allow wind tendencies all the way to the model bottom.
  tend_level = pver

  ! No spectrum; phase speed is just 0.
  c = 0._r8

  where( m2src < orom2min ) 
     tlb = mxdis
     tlb_level = src_level
  endwhere


end subroutine gw_rdg_src


!==========================================================================

subroutine gw_rdg_belowpeak(ncol, band, rdg_cd_llb, &
     t, mxdis, anixy, kwvrdg, zi, nm, ni, rhoi, &
     src_level , tau,  & 
     ubmsrc, nsrc, rsrc, m2src,tlb,bwv,Fr1,Fr2,Frx, & 
     tauoro,taudsw, hdspwv,hdspdw  )

  use gw_common, only: GWBand
  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Band to emit orographic waves in.
  ! Regardless, we will only ever emit into l = 0.
  type(GWBand), intent(in) :: band
  ! Drag coefficient for low-level flow
  real(r8), intent(in) :: rdg_cd_llb


  ! Midpoint temperatures. (K)
  real(r8), intent(in) :: t(ncol,pver)
  ! Height estimate for ridge (m) [anisotropic orography].
  real(r8), intent(in) :: mxdis(ncol)
  ! Anisotropy parameter [0-1] [anisotropic orography].
  real(r8), intent(in) :: anixy(ncol)
  ! Inverse cross-ridge lengthscale (m-1) [anisotropic orography].
  real(r8), intent(inout) :: kwvrdg(ncol)
  ! Interface altitudes above ground (m).
  real(r8), intent(in) :: zi(ncol,pver+1)
  ! Midpoint Brunt-Vaisalla frequencies (s-1).
  real(r8), intent(in) :: nm(ncol,pver)
  ! Interface Brunt-Vaisalla frequencies (s-1).
  real(r8), intent(in) :: ni(ncol,pver+1)
  ! Interface density (kg m-3).
  real(r8), intent(in) :: rhoi(ncol,pver+1)

  ! Indices of top gravity wave source level
  integer, intent(inout) :: src_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Top of low-level flow layer.
  real(r8), intent(inout) :: tlb(ncol)
  ! Bottom of linear wave region.
  real(r8), intent(inout) :: bwv(ncol)
  ! surface stress from linear waves.
  real(r8), intent(out) :: tauoro(ncol)
  ! surface stress for downslope wind regime.
  real(r8), intent(out) :: taudsw(ncol)

  ! Surface streamline displacement height for linear waves.
  real(r8), intent(out) :: hdspwv(ncol)
  ! Surface streamline displacement height for downslope wind regime.
  real(r8), intent(out) :: hdspdw(ncol)



  ! Froude numbers for flow/drag regimes
  real(r8), intent(in) :: Fr1(ncol), Fr2(ncol),Frx(ncol)

  ! Averages over source region.
  real(r8), intent(in) :: m2src(ncol) ! normalized non-hydro wavenumber
  real(r8), intent(in) :: nsrc(ncol)  ! B-V frequency.
  real(r8), intent(in) :: rsrc(ncol)  ! Density.
  real(r8), intent(in) :: ubmsrc(ncol) ! On-ridge wind.


  !logical, intent(in), optional :: forcetlb

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  real(r8) :: Coeff_LB(ncol),tausat,ubsrcx(ncol),dswamp
  real(r8) :: taulin(ncol),BetaMax

  ! ubsrcx introduced to account for situations with high shear, strong strat.
  do i = 1, ncol
        ubsrcx(i)    = max( ubmsrc(i)  , 0._r8 )
  end do

  do i = 1, ncol
     if ( m2src(i) > orom2min )   then 
        hdspwv(i) = min( mxdis(i) , Fr1(i) * ubsrcx(i) / nsrc(i) )
     else
        hdspwv(i) = 0._r8
     end if
  end do
  
  if (do_divstream) then
     do i = 1, ncol
        if ( m2src(i) > orom2min )   then 
           hdspdw(i) = min( mxdis(i) , Fr2(i) * ubsrcx(i) / nsrc(i) )
        else
           hdspdw(i) = 0._r8
        end if
     end do
  else
     do i = 1, ncol
        ! Needed only to mark where a DSW occurs
        if ( m2src(i) > orom2min )   then 
           hdspdw(i) = mxdis(i) 
        else
           hdspdw(i) = 0._r8
        end if
     end do
  end if

  ! Calculate form drag coefficient ("CD")
  !--------------------------------------
  Coeff_LB = rdg_cd_llb*anixy

  ! Determine the orographic c=0 source term following McFarlane (1987).
  ! Set the source top interface index to pver, if the orographic term is
  ! zero.
  ! 
  ! This formula is basically from
  !
  !      tau(src) = rho * u' * w'
  ! where 
  !      u' ~ N*h'  and w' ~ U*h'/b  (b="breite")
  !
  ! and 1/b has been replaced with k (kwvrdg) 
  !
  do i = 1, ncol
     if ( ( src_level(i) > 0 ) .and. ( m2src(i) > orom2min ) ) then
        tauoro(i) = kwvrdg(i) * ( hdspwv(i)**2 ) * rsrc(i) * nsrc(i) &
             * ubsrcx(i)
        taudsw(i) = kwvrdg(i) * ( hdspdw(i)**2 ) * rsrc(i) * nsrc(i) &
             * ubsrcx(i)
     else
        tauoro(i) = 0._r8
        taudsw(i) = 0._r8
     end if
  end do

  if (do_divstream) then
     do i = 1, ncol
           taulin(i) = 0._r8
     end do
  !---------------------------------------
  ! Need linear drag when divstream is not used is used
  !---------------------------------------
  else
     do i = 1, ncol
        if ( ( src_level(i) > 0 ) .and. ( m2src(i) > orom2min ) ) then
           taulin(i) = kwvrdg(i) * ( mxdis(i)**2 ) * rsrc(i) * nsrc(i) &
                * ubsrcx(i)
        else
           taulin(i) = 0._r8
        end if
     end do
  end if

  if ( do_divstream ) then
  ! Amplify DSW between Frx=1. and Frx=Frx1
     do i = 1,ncol
        dswamp=0._r8
        BetaMax   = C_BetaMax_DS * anixy(i)      
        if ( (Frx(i)>1._r8).and.(Frx(i)<=Frx1)) then
           dswamp = (Frx(i)-1._r8)*(Frx1-Frx(i))/(0.25_r8*(Frx1-1._r8)**2)
        end if
        taudsw(i) = (1._r8 + BetaMax*dswamp)*taudsw(i)
     end do
  else
  !-------------------
  ! Scinocca&McFarlane
  !--------------------
     do i = 1, ncol
        BetaMax   = C_BetaMax_SM * anixy(i)      
        if ( (Frx(i) >=1._r8) .and. (Frx(i) < 1.5_r8) ) then
           dswamp = 2._r8 * BetaMax * (Frx(i) -1._r8)
        else if ( ( Frx(i) >= 1.5_r8 ) .and. (Frx(i) < 3._r8 ) ) then
           dswamp = ( 1._r8 + BetaMax - (0.666_r8**2) ) * ( 0.666_r8*(3._r8 - Frx(i) ))**2  & 
                      + ( 1._r8 / Frx(i) )**2  -1._r8
        else
           dswamp    = 0._r8      
        end if
        if ( (Frx(i) >=1._r8) .and. (Frx(i) < 3._r8) ) then
          taudsw(i) = (1._r8 + dswamp )*taulin(i) - tauoro(i)
        else
          taudsw(i) = 0._r8   
        endif
        ! This code defines "taudsw" as SUM of freely-propagating
        ! DSW enhancement. Different than in SM2000
        taudsw(i) = taudsw(i) + tauoro(i) 
     end do
 !----------------------------------------------------
  end if

  
  do i = 1, ncol
     if ( m2src(i) > orom2min )   then 
        where ( ( zi(i,:) < mxdis(i) ) .and. ( zi(i,:) >= bwv(i) ) )
             tau(i,0,:) =  tauoro(i)
        else where ( ( zi(i,:) < bwv(i) ) .and. ( zi(i,:) >= tlb(i) ) )
             tau(i,0,:) =  tauoro(i) +( taudsw(i)-tauoro(i) )* &
                                         ( bwv(i) - zi(i,:) ) / &
                                         ( bwv(i) - tlb(i) )
        endwhere
        ! low-level form drag on obstacle. Quantity kwvrdg (~1/b) appears for consistency
        ! with tauoro and taudsw forms. Should be weighted by L*b/A_g before applied to flow.
        where ( ( zi(i,:) < tlb(i) ) .and. ( zi(i,:) >= 0._r8 ) )
             tau(i,0,:) =  taudsw(i) +  &
                           Coeff_LB(i) * kwvrdg(i) * rsrc(i) * 0.5_r8 * (ubsrcx(i)**2) * ( tlb(i) - zi(i,:) )
        endwhere
 
        if (do_smooth_regimes) then
        !  This blocks accounts for case where both mxdis and tlb fall
        !  between adjacent edges
           do k=1,pver
              if ( (zi(i,k) >= tlb(i)).and.(zi(i,k+1) < tlb(i)).and. &
                   (zi(i,k) >= mxdis(i)).and.(zi(i,k+1) < mxdis(i)) ) then
                 src_level(i) = src_level(i)-1
                 tau(i,0,k) = tauoro(i)
              end if
           end do
        end if 

     else     !----------------------------------------------
             ! This block allows low-level dynamics to occur
             ! even if m2 is less than orom2min
        where ( ( zi(i,:) < tlb(i) ) .and. ( zi(i,:) >= 0._r8 ) )
               tau(i,0,:) =  taudsw(i) +  &
                   Coeff_LB(i) * kwvrdg(i) * rsrc(i) * 0.5_r8 * &
                   (ubsrcx(i)**2) * ( tlb(i) - zi(i,:) )
        endwhere
     endif
  end do

  ! This may be redundant with newest version of gw_drag_prof.
  ! That code reaches down to level k=src_level+1. (jtb 1/5/16)
  do i = 1, ncol
     k=src_level(i)
     if ( ni(i,k) > orostratmin ) then
         tausat    =  (Fr_c**2) * kwvrdg(i) * rhoi(i,k) * ubsrcx(i)**3 / &
              (1._r8*ni(i,k)) 
     else
         tausat = 0._r8
     endif 
     tau(i,0,src_level(i)) = min( tauoro(i), tausat ) 
  end do



  ! Final clean-up. Do nothing if obstacle less than orohmin
  do i = 1, ncol
     if ( mxdis(i) < orohmin ) then
        tau(i,0,:) = 0._r8 
        tauoro(i)  = 0._r8
        taudsw(i)  = 0._r8
     endif 
  end do

          ! Disable vertical propagation if Scorer param is 
          ! too small.
  do i = 1, ncol
     if ( m2src(i) <= orom2min ) then
        src_level(i)=1
     endif 
  end do



end subroutine gw_rdg_belowpeak

!==========================================================================
subroutine gw_rdg_break_trap(ncol, band, &
     zi, nm, ni, ubm, ubi, rhoi, kwvrdg, bwv, tlb, wbr, & 
     src_level, tlb_level, & 
     hdspwv, hdspdw, mxdis, &
     tauoro, taudsw,  tau, & 
     ldo_trapped_waves, wdth_kwv_scale_in )
  use gw_common, only: GWBand
  !-----------------------------------------------------------------------
  ! Parameterization of high-drag regimes and trapped lee-waves for CAM
  !
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Band to emit orographic waves in.
  ! Regardless, we will only ever emit into l = 0.
  type(GWBand), intent(in) :: band


  ! Height estimate for ridge (m) [anisotropic orography].
  !real(r8), intent(in) :: mxdis(ncol)
  ! Horz wavenumber for ridge (1/m) [anisotropic orography].
  real(r8), intent(in) :: kwvrdg(ncol)
  ! Interface altitudes above ground (m).
  real(r8), intent(in) :: zi(ncol,pver+1)
  ! Midpoint Brunt-Vaisalla frequencies (s-1).
  real(r8), intent(in) :: nm(ncol,pver)
  ! Interface Brunt-Vaisalla frequencies (s-1).
  real(r8), intent(in) :: ni(ncol,pver+1)

  ! Indices of gravity wave sources.
  integer, intent(inout) :: src_level(ncol), tlb_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Wave Reynolds stresses at source.
  real(r8), intent(inout) :: taudsw(ncol),tauoro(ncol)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(in) :: ubm(ncol,pver)
  real(r8), intent(in) :: ubi(ncol,pver+1)
  ! Interface density (kg m-3).
  real(r8), intent(in) :: rhoi(ncol,pver+1)

  ! Top of low-level flow layer.
  real(r8), intent(in) :: tlb(ncol)
  ! Bottom of linear wave region.
  real(r8), intent(in) :: bwv(ncol)

  ! Surface streamline displacement height for linear waves.
  real(r8), intent(in) :: hdspwv(ncol)
  ! Surface streamline displacement height for downslope wind regime.
  real(r8), intent(in) :: hdspdw(ncol)
  ! Ridge height.
  real(r8), intent(in) :: mxdis(ncol)


  ! Wave breaking level
  real(r8), intent(out) :: wbr(ncol)

  logical, intent(in), optional :: ldo_trapped_waves
  real(r8), intent(in), optional :: wdth_kwv_scale_in

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k, kp1, non_hydro
  real(r8):: m2(ncol,pver),delz(ncol),tausat(ncol),trn(ncol)
  real(r8):: wbrx(ncol)
  real(r8):: phswkb(ncol,pver+1)
  logical :: lldo_trapped_waves
  real(r8):: wdth_kwv_scale
  ! Indices of important levels.
  integer :: trn_level(ncol)

  if (present(ldo_trapped_waves)) then
     lldo_trapped_waves = ldo_trapped_waves
     if(lldo_trapped_waves) then
       non_hydro = 1
     else
       non_hydro = 0
     endif
  else
     lldo_trapped_waves = .false.
     non_hydro = 0
  endif

  if (present(wdth_kwv_scale_in)) then
     wdth_kwv_scale = wdth_kwv_scale_in
  else
     wdth_kwv_scale = 1._r8
  endif

  ! Calculate vertical wavenumber**2
  !---------------------------------
  m2 = (nm  / (abs(ubm)+.01_r8))**2
  do k=pver,1,-1
     m2(:,k) = m2(:,k) - non_hydro*(wdth_kwv_scale*kwvrdg)**2
     ! sweeping up, zero out m2 above first occurence
     ! of m2(:,k)<=0
     kp1=min( k+1, pver )
     where( (m2(:,k) <= 0.0_r8 ).or.(m2(:,kp1) <= 0.0_r8 ) )
        m2(:,k) = 0._r8
     endwhere
  end do

  ! Take square root of m**2 and 
  ! do vertical integral to find
  ! WKB phase.
  !-----------------------------
  m2 = SQRT( m2 )
  phswkb(:,:)=0
  do k=pver,1,-1
     where( zi(:,k) > tlb(:) )
        delz(:) = min( zi(:,k)-zi(:,k+1) , zi(:,k)-tlb(:) ) 
        phswkb(:,k) = phswkb(:,k+1) + m2(:,k)*delz(:) 
     endwhere
  end do

  ! Identify top edge of layer in which phswkb reaches 3*pi/2
  ! - approximately the "breaking level"
  !----------------------------------------------------------
  wbr(:)=0._r8
  wbrx(:)=0._r8
  if (do_smooth_regimes) then
     do k=pver,1,-1
     where( (phswkb(:,k+1)<1.5_r8*pii).and.(phswkb(:,k)>=1.5_r8*pii) & 
            .and.(hdspdw(:)>hdspwv(:)) )
        wbr(:)  = zi(:,k)  
        ! Extrapolation to make regime
        ! transitions smoother
        wbrx(:) = zi(:,k)   - ( phswkb(:,k) -  1.5_r8*pii ) &
                            / ( m2(:,k) + 1.e-6_r8 )
        src_level(:) = k-1
     endwhere
     end do
  else
     do k=pver,1,-1
     where( (phswkb(:,k+1)<1.5_r8*pii).and.(phswkb(:,k)>=1.5_r8*pii) & 
            .and.(hdspdw(:)>hdspwv(:)) )
        wbr(:)  = zi(:,k)
        src_level(:) = k
     endwhere
     end do
  end if

  ! Adjust tauoro at new source levels if needed.
  ! This is problematic if Fr_c<1.0. Not sure why.
  !----------------------------------------------------------
  if (do_adjust_tauoro) then 
     do i = 1,ncol
        if (wbr(i) > 0._r8 ) then
            tausat(i) = (Fr_c**2) * kwvrdg(i)  * rhoi( i, src_level(i) ) & 
                      * abs(ubi(i , src_level(i) ))**3  &
                      / ni( i , src_level(i) ) 
            tauoro(i) = min( tauoro(i), tausat(i) )
        end if
     end do
  end if

  if (do_smooth_regimes) then
     do i = 1, ncol
     do k=1,pver+1
        if ( ( zi(i,k) <= wbr(i) ) .and. ( zi(i,k) > tlb(i) ) ) then
           tau(i,0,k) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                          ( wbrx(i) - zi(i,k) ) / &
                          ( wbrx(i) - tlb(i)  )
           tau(i,0,k) = max( tau(i,0,k), tauoro(i) ) 
        endif
     end do   
     end do
  else
  ! Following is for backwards B4B compatibility with earlier versions
  ! ("N1" and "N5" -- Note: "N5" used do_backward_compat=.true.)
     if (.not.do_backward_compat) then
        do i = 1, ncol
        do k=1,pver+1
           if ( ( zi(i,k) <  wbr(i) ) .and. ( zi(i,k) >= tlb(i) ) ) then
              tau(i,0,k) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                            ( wbr(i) - zi(i,k) ) / &
                            ( wbr(i) - tlb(i)  )
           endif
        end do   
        end do
     else
        do i = 1, ncol
        do k=1,pver+1
           if ( ( zi(i,k) <= wbr(i) ) .and. ( zi(i,k) > tlb(i) ) ) then
              tau(i,0,k) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                            ( wbr(i) - zi(i,k) ) / &
                            ( wbr(i) - tlb(i)  )
           endif
        end do   
        end do
     end if
  end if
  
  if (lldo_trapped_waves) then 
     
  ! Identify top edge of layer in which Scorer param drops below 0
  ! - approximately the "turning level"
  !----------------------------------------------------------
     trn(:)=1.e8_r8
     trn_level(:) = 0 ! pver+1
     where( m2(:,pver)<= 0._r8 )
         trn(:) = zi(:,pver)
         trn_level(:) = pver
     endwhere
     do k=pver-1,1,-1
        where( (m2(:,k+1)> 0._r8).and.(m2(:,k)<= 0._r8) )
           trn(:) = zi(:,k)
           trn_level(:) = k
        endwhere
     end do

     do i = 1,ncol
     ! Case: Turning below mountain top
        if ( (trn(i) < mxdis(i)).and.(trn_level(i)>=1) ) then
            tau(i,0,:) =  tau(i,0,:) - max( tauoro(i),taudsw(i) )
            tau(i,0,:) =  max( tau(i,0,:) , 0._r8 )
            tau(i,0,1:tlb_level(i))=0._r8
            src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning but no breaking
        if ( (wbr(i) == 0._r8 ).and.(trn(i)>mxdis(i)).and.(trn_level(i)>=1) ) then
           where ( ( zi(i,:) <= trn(i) ) .and. ( zi(i,:) >= bwv(i) ) )
               tau(i,0,:) =  tauoro(i) * &
                             ( trn(i) - zi(i,:) ) / &
                             ( trn(i) - bwv(i)  )
           end where
           src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning AND breaking. Turning ABOVE breaking
        if ( (wbr(i) > 0._r8 ).and.(trn(i) >= wbr(i)).and.(trn_level(i)>=1) ) then
           where ( ( zi(i,:) <= trn(i) ) .and. ( zi(i,:) >= wbr(i) ) )
               tau(i,0,:) =   tauoro(i) * &
                             ( trn(i) - zi(i,:) ) / &
                             ( trn(i) - wbr(i)  )
           endwhere
           src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning AND breaking. Turning BELOW breaking
        if ( (wbr(i) > 0._r8 ).and.(trn(i) < wbr(i)).and.(trn_level(i)>=1) ) then
           tauoro(i) = 0._r8
           where ( ( zi(i,:) < wbr(i) ) .and. ( zi(i,:) >= tlb(i) ) )
               tau(i,0,:) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                             ( wbr(i) - zi(i,:) ) / &
                             ( wbr(i) - tlb(i)  )
           endwhere
           src_level(i) = 1 ! disable any more tau calculation
        end if
     end do
  end if

  end subroutine gw_rdg_break_trap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module gw_rdg
