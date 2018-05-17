module diag_module

! diagnostic calcs
!
! REVISION HISTORY:
!   05.09.10  Rasch   Creation of compute_vdot_gradp
!   05.10.18  Sawyer  Revisions for 2D decomp, placed in module
!   07.01.29  Chen    Removed pft2d calculation for OMGA (is in cd_core)

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use pmgrid,         only: plon, plev, plevp
use cam_logfile,    only: iulog
use cam_history,    only: addfld, outfld, add_default, horiz_only
use cam_abortutils, only: endrun

implicit none
private
save

public :: &
   fv_diag_init,       &
   fv_diag_am_calc,    &
   compute_vdot_gradp

real(r8) :: rplon
real(r8) :: iref_p(plevp)              ! interface reference pressure for vertical interpolation
integer  :: ip_b                       ! level index where hybrid levels become purely pressure
integer  :: zm_limit

!========================================================================================
CONTAINS
!========================================================================================

subroutine fv_diag_init()

   use hycoef,       only : hyai, hybi, ps0

   ! local variables
   integer :: k
   logical :: history_waccm
   !---------------------------------------------------------------------------

   rplon    = 1._r8/plon
   zm_limit = plon/3

   !-------------------------------------------------------------
   ! Calculate reference pressure
   !-------------------------------------------------------------
   do k = 1, plevp
      iref_p(k) = (hyai(k) + hybi(k)) * ps0
   end do
   if( masterproc ) then
      write(iulog,*) 'fv_diag_inti: iref_p'
      write(iulog,'(1p5g15.7)') iref_p(:)
   end if

   !-------------------------------------------------------------
   ! Find level where hybrid levels become purely pressure 
   !-------------------------------------------------------------
   ip_b = -1
   do k = 1,plev
      if( hybi(k) == 0._r8 ) ip_b = k
   end do

   ! Fields for diagnosing angular momentum conservation.  They are supplemental
   ! to the fields computed by do_circulation_diags
   
   call addfld ('dUzm' ,(/ 'ilev' /),'A','M/S','Zonal-Mean U dyn increm - defined on ilev', gridname='fv_centers_zonal' )
   call addfld ('dVzm' ,(/ 'ilev' /),'A','M/S','Zonal-Mean V dyn increm - defined on ilev', gridname='fv_centers_zonal' )
   call addfld ('dUazm',(/ 'ilev' /),'A','M/S','Zonal-Mean U adv increm - defined on ilev', gridname='fv_centers_zonal' )
   call addfld ('dVazm',(/ 'ilev' /),'A','M/S','Zonal-Mean V adv increm - defined on ilev', gridname='fv_centers_zonal' )
   call addfld ('dUfzm',(/ 'ilev' /),'A','M/S','Zonal-Mean U fixer incr - defined on ilev', gridname='fv_centers_zonal' )
   call addfld ('dU',   (/ 'lev' /), 'A','K',  'U dyn increm', gridname='fv_centers' )
   call addfld ('dV',   (/ 'lev' /), 'A','K',  'V dyn increm', gridname='fv_centers' )
   call addfld ('dUa',  (/ 'lev' /), 'A','K',  'U adv increm', gridname='fv_centers' )
   call addfld ('dVa',  (/ 'lev' /), 'A','K',  'V adv increm', gridname='fv_centers' )
   call addfld ('dUf',  (/ 'lev' /), 'A','K',  'U fixer incr', gridname='fv_centers' )
    
   call add_default ('dUzm'  ,1, ' ')
   call add_default ('dVzm'  ,1, ' ')
   call add_default ('dUazm' ,1, ' ')
   call add_default ('dVazm' ,1, ' ')
   call add_default ('dUfzm' ,1, ' ')
   call add_default ('dU' ,   1, ' ')
   call add_default ('dV' ,   1, ' ')
   call add_default ('dUa',   1, ' ')
   call add_default ('dVa',   1, ' ')
   call add_default ('dUf',   1, ' ')

end subroutine fv_diag_init

!========================================================================================

subroutine fv_diag_am_calc(grid, ps, pe, du3, dv3, dua3, dva3, duf3)

   ! Compute fields for diagnosing angular momentum conservation.  They are supplemental
   ! to the fields computed by do_circulation_diags

    use dynamics_vars, only          : T_FVDYCORE_GRID
    use interpolate_data, only       : vertinterp
#ifdef SPMD
    use mpishorthand,       only     : mpilog, mpiint
    use parutilitiesmodule, only     : pargatherint
#endif

!-------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------
    type(T_FVDYCORE_GRID), intent(in) :: grid                        ! FV Dynamics grid
    real(r8), intent(in)  :: ps(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)         ! surface pressure (pa)
    real(r8), intent(in)  :: pe(grid%ifirstxy:grid%ilastxy,plevp,grid%jfirstxy:grid%jlastxy)   ! interface pressure (pa)
    real(r8), intent(in)  :: du3 (grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)  ! U increment, total (m/s/timestep)
    real(r8), intent(in)  :: dv3 (grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)  ! V increment, total (m/s/timestep)
    real(r8), intent(in)  :: dua3(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)  ! U increment, advec (m/s/timestep)
    real(r8), intent(in)  :: dva3(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)  ! V increment, advec (m/s/timestep)
    real(r8), intent(in)  :: duf3(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)  ! U increment, fixer (m/s/timestep)

!-------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------
    
    real(r8) :: pinterp
    real(r8) :: pm(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! mid-point pressure
    real(r8) :: psurf

    real(r8) :: dui (grid%ifirstxy:grid%ilastxy,plevp)       ! interp. zonal mean U increment, total FV 
    real(r8) :: dvi (grid%ifirstxy:grid%ilastxy,plevp)       ! interp. zonal mean V increment, total FV 
    real(r8) :: duai(grid%ifirstxy:grid%ilastxy,plevp)       ! interp. zonal mean U increment, advection
    real(r8) :: dvai(grid%ifirstxy:grid%ilastxy,plevp)       ! interp. zonal mean V increment, advection
    real(r8) :: dufi(grid%ifirstxy:grid%ilastxy,plevp)       ! interp. zonal mean U increment, fixer    
    
    real(r8) :: dum (plevp)                                  ! zonal mean U increment, total FV 
    real(r8) :: dvm (plevp)                                  ! zonal mean V increment, total FV 
    real(r8) :: duam(plevp)                                  ! zonal mean U increment, advection
    real(r8) :: dvam(plevp)                                  ! zonal mean V increment, advection
    real(r8) :: dufm(plevp)                                  ! zonal mean U increment, fixer    

    real(r8) :: rdiv(plevp)
    
    integer  :: ip_gm1g(plon,grid%jfirstxy:grid%jlastxy)     ! contains level index-1 where blocked points begin
    integer  :: zm_cnt(plevp)                                ! counter
    integer  :: i,j,k
    integer  :: nlons

    logical  :: has_zm(plevp,grid%jfirstxy:grid%jlastxy)                   ! .true. the (z,y) point is a valid zonal mean 
    integer  :: ip_gm1(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy) ! contains level index-1 where blocked points begin

    real(r8) :: du2d (plevp,grid%jfirstxy:grid%jlastxy)              ! zonally averaged U dycore increment
    real(r8) :: dv2d (plevp,grid%jfirstxy:grid%jlastxy)              ! zonally averaged V dycore increment
    real(r8) :: dua2d(plevp,grid%jfirstxy:grid%jlastxy)              ! zonally averaged U advect increment
    real(r8) :: dva2d(plevp,grid%jfirstxy:grid%jlastxy)              ! zonally averaged V advect increment
    real(r8) :: duf2d(plevp,grid%jfirstxy:grid%jlastxy)              ! zonally averaged U fixer  increment

    real(r8) :: tmp2(grid%ifirstxy:grid%ilastxy)  
    real(r8) :: tmp3(grid%ifirstxy:grid%ilastxy,plevp)  
    real(r8) :: tmph(grid%ifirstxy:grid%ilastxy,plev)  

    integer :: beglat, endlat                          ! begin,end latitude indicies
    integer :: beglon, endlon                          ! begin,end longitude indicies

    beglon = grid%ifirstxy
    endlon = grid%ilastxy
    beglat = grid%jfirstxy
    endlat = grid%jlastxy

!omp parallel do private (i,j,k,psurf)
lat_loop1 : &
    do j = beglat, endlat
       do k = 1, plev
          do i = beglon, endlon
             pm(i,k,j) = 0.5_r8 * ( pe(i,k,j) + pe(i,k+1,j) )
          end do
       end do
!-------------------------------------------------------------
! Keep track of where the bottom is in each column 
! (i.e., largest index for which P(k) <= PS)
!-------------------------------------------------------------
       ip_gm1(:,j) = plevp 
       do i = beglon, endlon
          psurf = ps(i,j)
          do k = ip_b+1, plevp
             if( iref_p(k) <= psurf ) then
                ip_gm1(i,j) = k
             end if
          end do
       end do
    end do lat_loop1

    nlons = endlon - beglon + 1

#ifdef SPMD    
    if( grid%twod_decomp == 1 ) then
       if (grid%iam .lt. grid%npes_xy) then
          call pargatherint( grid%commxy_x, 0, ip_gm1, grid%strip2dx, ip_gm1g )
       endif
    else
       ip_gm1g(:,:) = ip_gm1(:,:)
    end if
#else
    ip_gm1g(:,:) = ip_gm1(:,:)
#endif
#ifdef SPMD    
    if( grid%myidxy_x == 0 ) then
#endif
lat_loop2 : &
       do j = beglat, endlat
          zm_cnt(:ip_b) = plon
          do k = ip_b+1, plevp
             zm_cnt(k) = count( ip_gm1g(:,j) >= k )
          end do
          has_zm(:ip_b,j) = .true.
          do k = ip_b+1, plevp
             has_zm(k,j) = zm_cnt(k) >= zm_limit
          end do
       end do lat_loop2
#ifdef SPMD    
    end if
    if( grid%twod_decomp == 1 ) then
       call mpibcast( has_zm, plevp*(endlat-beglat+1), mpilog, 0, grid%commxy_x )
       call mpibcast( ip_gm1g, plon*(endlat-beglat+1), mpiint, 0, grid%commxy_x )
    end if
#endif

lat_loop3 : &
    do j = beglat, endlat
!-------------------------------------------------------------
! Vertical interpolation
!-------------------------------------------------------------
       do k = 1, plevp
          pinterp = iref_p(k)
!-------------------------------------------------------------
!      Zonal & meridional velocity increments
!-------------------------------------------------------------
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           du3(beglon,1,j) , dui (beglon,k) )
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           dv3(beglon,1,j) , dvi (beglon,k) )
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           dua3(beglon,1,j), duai(beglon,k) )
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           dva3(beglon,1,j), dvai(beglon,k) )
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           duf3(beglon,1,j), dufi(beglon,k) )
       end do

!-------------------------------------------------------------
! Calculate zonal averages
!-------------------------------------------------------------
       do k = ip_b+1, plevp
          where( ip_gm1(beglon:endlon,j) < k )
             dui (beglon:endlon,k)= 0._r8
             dvi (beglon:endlon,k)= 0._r8
             duai(beglon:endlon,k)= 0._r8
             dvai(beglon:endlon,k)= 0._r8
             dufi(beglon:endlon,k)= 0._r8
          endwhere
       end do

       call par_xsum(grid, dui,  plevp, dum)
       call par_xsum(grid, dvi,  plevp, dvm)
       call par_xsum(grid, duai, plevp, duam)
       call par_xsum(grid, dvai, plevp, dvam)
       call par_xsum(grid, dufi, plevp, dufm)

       do k = 1, ip_b
          du2d(k,j) = dum(k) * rplon
          dv2d(k,j) = dvm(k) * rplon
          dua2d(k,j)= duam(k)* rplon
          dva2d(k,j)= dvam(k)* rplon
          duf2d(k,j)= dufm(k)* rplon
       end do

       do k = ip_b+1, plevp
          if( has_zm(k,j) ) then
             rdiv(k)   = 1._r8/count( ip_gm1g(:,j) >= k )
             du2d(k,j) = dum(k) * rdiv(k)
             dv2d(k,j) = dvm(k) * rdiv(k)
             dua2d(k,j)= duam(k)* rdiv(k)
             dva2d(k,j)= dvam(k)* rdiv(k)
             duf2d(k,j)= dufm(k)* rdiv(k)
          else
             du2d(k,j) = 0._r8
             dv2d(k,j) = 0._r8
             dua2d(k,j)= 0._r8
             dva2d(k,j)= 0._r8
             duf2d(k,j)= 0._r8
          end if
       end do

    end do lat_loop3

!-------------------------------------------------------------
! Do the output
!-------------------------------------------------------------
    latloop: do j = beglat,endlat

!-------------------------------------------------------------
! zonal-mean output
!-------------------------------------------------------------

       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = du2d(k,j)
       enddo
       call outfld( 'dUzm', tmp3(grid%ifirstxy,:), 1, j )

       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = dv2d(k,j)
       enddo
       call outfld( 'dVzm', tmp3(grid%ifirstxy,:), 1, j )

       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = dua2d(k,j)
       enddo
       call outfld( 'dUazm', tmp3(grid%ifirstxy,:), 1, j )

       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = dva2d(k,j)
       enddo
       call outfld( 'dVazm', tmp3(grid%ifirstxy,:), 1, j )

       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = duf2d(k,j)
       enddo
       call outfld( 'dUfzm', tmp3(grid%ifirstxy,:), 1, j )

!-------------------------------------------------------------
! 3D output
!-------------------------------------------------------------

       do k = 1,plev
          do i = beglon,endlon
             tmph(i,k) = du3(i,k,j)
          enddo
       enddo
       call outfld( 'dU', tmph, nlons, j )

       do k = 1,plev
          do i = beglon,endlon
             tmph(i,k) = dv3(i,k,j)
          enddo
       enddo
       call outfld( 'dV', tmph, nlons, j )

       do k = 1,plev
          do i = beglon,endlon
             tmph(i,k) = dua3(i,k,j)
          enddo
       enddo
       call outfld( 'dUa', tmph, nlons, j )

       do k = 1,plev
          do i = beglon,endlon
             tmph(i,k) = dva3(i,k,j)
          enddo
       enddo
       call outfld( 'dVa', tmph, nlons, j )

       do k = 1,plev
          do i = beglon,endlon
             tmph(i,k) = duf3(i,k,j)
          enddo
       enddo
       call outfld( 'dUf', tmph, nlons, j )

    enddo latloop

end subroutine fv_diag_am_calc

!========================================================================================

  subroutine compute_vdot_gradp(grid, dt, frac, cx, cy, pexy, omgaxy )

  use shr_kind_mod, only :  r8 => shr_kind_r8
  use dynamics_vars, only : T_FVDYCORE_GRID
#if defined( SPMD )
  use mod_comm, only: mp_send3d, mp_recv3d, &
                      mp_sendirr, mp_recvirr
#endif

  implicit none

! !INPUT PARAMETERS:
  type (T_FVDYCORE_GRID), intent(in) :: grid
  real(r8), intent(in):: dt
  real(r8), intent(in):: frac

  real(r8), intent(in):: cx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
  real(r8), intent(in):: cy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)
  real(r8), target, intent(in)::   &
    pexy(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy) ! P (pascal) at layer edges
  real(r8), target, intent(inout):: &
    omgaxy(grid%ifirstxy:grid%ilastxy,grid%km,grid%jfirstxy:grid%jlastxy) ! vert. press. velocity (pa/sec)

! Local 
  integer :: im       ! dimension in east-west
  integer :: jm       ! dimension in North-South
  integer :: km       ! number of Lagrangian layers
  integer :: jfirst   ! starting latitude index for MPI
  integer :: jlast    ! ending latitude index for MPI
  integer :: kfirst   ! starting level index for MPI
  integer :: klast    ! ending level index for MPI
  integer :: js2g0    ! j==1 not included
  integer :: jn2g0    ! j==jm not included

  real(r8) :: pm(grid%im, grid%jfirst-1:grid%jlast+1)
  real(r8) :: grad(grid%im, grid%jfirst:grid%jlast+1)
  real(r8) :: fac, sum1

  real(r8), pointer :: pe(:,:,:)      ! YZ version of edge pressures
  real(r8), pointer :: omga(:,:,:)    ! YZ version of vert. vel.

  real(r8), parameter :: half =  0.5_r8
  real(r8), parameter :: zero =  0.0_r8

  integer  :: i,j,k

#if defined( SPMD )
  integer  :: iam, dest, src, npr_y, npes_yz
  real(r8) :: penorth(grid%im, grid%kfirst:grid%klast+1)
  real(r8) :: pesouth(grid%im, grid%kfirst:grid%klast+1)
#endif

  im     = grid%im
  jm     = grid%jm
  km     = grid%km
  jfirst = grid%jfirst
  jlast  = grid%jlast
  kfirst = grid%kfirst
  klast  = grid%klast
  js2g0  = grid%js2g0
  jn2g0  = grid%jn2g0

  fac = half / (dt * frac)

#if defined( SPMD )
  if ( grid%twod_decomp == 1 ) then
    allocate(pe(im,kfirst:klast+1,jfirst:jlast))
    allocate(omga(im,kfirst:klast,jfirst:jlast))
    call mp_sendirr( grid%commxy, grid%ikj_xy_to_yz%SendDesc,                  &
                     grid%ikj_xy_to_yz%RecvDesc, omgaxy, omga,                 &
                     modc=grid%modc_dynrun )
    call mp_recvirr( grid%commxy, grid%ikj_xy_to_yz%SendDesc,                  &
                     grid%ikj_xy_to_yz%RecvDesc, omgaxy, omga,                 &
                     modc=grid%modc_dynrun )
    call mp_sendirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                     grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                     modc=grid%modc_dynrun )
    call mp_recvirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                     grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                     modc=grid%modc_dynrun )
  else
    pe => pexy
    omga => omgaxy
  endif
  iam   = grid%iam
  npes_yz   = grid%npes_yz
 if (iam .lt. npes_yz) then
  npr_y = grid%npr_y
  dest  = iam+1
  src   = iam-1
  if ( mod(iam+1,npr_y) == 0 ) dest = -1
  if ( mod(iam,npr_y) == 0 ) src = -1

!
! Have to give more thought to the source and destination for 2D
!
  call mp_send3d(grid%commyz, dest, src, im, km+1, jm,       &
                 1, im, kfirst, klast+1, jfirst, jlast,     &
                 1, im, kfirst, klast+1, jlast, jlast, pe)
  call mp_recv3d(grid%commyz, src, im, km+1, jm,              &
                 1, im, kfirst, klast+1, jfirst-1, jfirst-1, &
                 1, im, kfirst, klast+1, jfirst-1, jfirst-1, pesouth)
  call mp_send3d(grid%commyz, src, dest, im, km+1, jm,       &
                 1, im, kfirst, klast+1, jfirst, jlast,     &
                 1, im, kfirst, klast+1, jfirst, jfirst, pe)
  call mp_recv3d(grid%commyz, dest, im, km+1, jm,            &
                 1, im, kfirst, klast+1, jlast+1, jlast+1,  &
                 1, im, kfirst, klast+1, jlast+1, jlast+1, penorth)
 end if  !  (iam .lt. npes_yz)
#else
  pe => pexy
  omga => omgaxy
#endif

!$omp parallel do private(i,j,k,pm,grad, sum1)
  do k=kfirst,klast

! Compute layer mean p
     do j=jfirst,jlast
        do i=1,im
           pm(i,j) = half * ( pe(i,k,j) + pe(i,k+1,j) )
        enddo
     enddo

#if defined( SPMD )
     if ( jfirst/=1 ) then
        do i=1,im
           pm(i,jfirst-1) = half * ( pesouth(i,k) + pesouth(i,k+1))
        enddo
     endif

     if ( jlast/=jm ) then
        do i=1,im
           pm(i,jlast+1) = half * ( penorth(i,k) + penorth(i,k+1))
        enddo
     endif
#endif

     do j=js2g0,jn2g0
           i=1
           grad(i,j) = fac * cx(i,j,k) * (pm(i,j)-pm(im,j)) 
        do i=2,im
           grad(i,j) = fac * cx(i,j,k) * (pm(i,j)-pm(i-1,j)) 
        enddo
     enddo

     do j=js2g0,jn2g0
        do i=1,im-1
           omga(i,k,j) = omga(i,k,j) + grad(i,j) + grad(i+1,j)
        enddo
           i=im
           omga(i,k,j) = omga(i,k,j) + grad(i,j) + grad(1,j)
     enddo

     do j=js2g0,min(jm,jlast+1)
        do i=1,im
           grad(i,j) = fac * cy(i,j,k) * (pm(i,j)-pm(i,j-1)) 
        enddo
     enddo

     do j=js2g0,jn2g0
        do i=1,im
           omga(i,k,j) = omga(i,k,j) + grad(i,j) + grad(i,j+1)
        enddo
     enddo

! Note: Since V*grad(P) at poles are harder to compute accurately we use the average of sourding points
!       to be used as input to physics.

     if ( jfirst==1 ) then
        sum1 = zero
        do i=1,im
           sum1 = sum1 + omga(i,k,2)
        enddo
        sum1 = sum1 / real(im,r8)
        do i=1,im
           omga(i,k,1) = sum1
        enddo
     endif

     if ( jlast==jm ) then
        sum1 = zero
        do i=1,im
           sum1 = sum1 + omga(i,k,jm-1)
        enddo
        sum1 = sum1 / real(im,r8)
        do i=1,im
           omga(i,k,jm) = sum1
        enddo
     endif
  enddo

#if defined( SPMD)
  if ( grid%twod_decomp == 1 ) then
!
! Transpose back to XY  (if 1D, the changes to omgaxy were made in place)
!
    call mp_sendirr( grid%commxy, grid%ikj_yz_to_xy%SendDesc,                  &
                     grid%ikj_yz_to_xy%RecvDesc, omga, omgaxy,                 &
                     modc=grid%modc_dynrun )
    call mp_recvirr( grid%commxy, grid%ikj_yz_to_xy%SendDesc,                  &
                     grid%ikj_yz_to_xy%RecvDesc, omga, omgaxy,                 &
                     modc=grid%modc_dynrun )
    deallocate( pe )
    deallocate( omga )
  endif
#endif

  end subroutine compute_vdot_gradp

end module diag_module
