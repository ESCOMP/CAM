!-----------------------------------------------------------------------------
! circulation diagnostics -- terms of the Transformed Eulerian Mean (TEM) equation
!-----------------------------------------------------------------------------
module ctem

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use pmgrid,       only: plon, plev, plevp
  use cam_logfile,  only: iulog
  use cam_history,  only: addfld, outfld, add_default, horiz_only
  use cam_abortutils, only: endrun

  implicit none
  private
  
  public :: ctem_readnl
  public :: ctem_init
  public :: ctem_diags
  public :: do_circulation_diags

  real(r8) :: rplon
  real(r8) :: iref_p(plevp)              ! interface reference pressure for vertical interpolation
  integer  :: ip_b                       ! level index where hybrid levels become purely pressure
  integer  :: zm_limit

  logical :: do_circulation_diags = .false.

contains

!================================================================================

  subroutine ctem_diags( u3, v3, omga, pt, h2o, ps, pe, grid)

    use physconst, only              : zvir, cappa
    use dynamics_vars, only          : T_FVDYCORE_GRID
    use hycoef, only                 : ps0
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
    real(r8), intent(in)  :: u3(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)    ! zonal velocity (m/s)
    real(r8), intent(in)  :: v3(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)    ! meridional velocity (m/s)
    real(r8), intent(in)  :: omga(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)  ! pressure velocity
    real(r8), intent(in)  :: pe(grid%ifirstxy:grid%ilastxy,plevp,grid%jfirstxy:grid%jlastxy)   ! interface pressure (pa)
    real(r8), intent(in)  :: pt(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,plev)    ! virtual temperature
    real(r8), intent(in)  :: h2o(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,plev)   ! water constituent (kg/kg)

!-------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------
    real(r8), parameter :: hscale = 7000._r8          ! pressure scale height
    real(r8), parameter :: navp   = 1.e35_r8
    
    real(r8) :: pinterp
    real(r8) :: w(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)          ! vertical velocity
    real(r8) :: th(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! pot. temperature

    real(r8) :: pm(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! mid-point pressure

    real(r8) :: pexf     ! Exner function
    real(r8) :: psurf

    real(r8) :: ui(grid%ifirstxy:grid%ilastxy,plevp)         ! interpolated zonal velocity
    real(r8) :: vi(grid%ifirstxy:grid%ilastxy,plevp)         ! interpolated meridional velocity
    real(r8) :: wi(grid%ifirstxy:grid%ilastxy,plevp)         ! interpolated vertical velocity
    real(r8) :: thi(grid%ifirstxy:grid%ilastxy,plevp)        ! interpolated pot. temperature
    
    real(r8) :: um(plevp)                                    ! zonal mean zonal velocity
    real(r8) :: vm(plevp)                                    ! zonal mean meridional velocity
    real(r8) :: wm(plevp)                                    ! zonal mean vertical velocity
    real(r8) :: thm(plevp)                                   ! zonal mean pot. temperature

    real(r8) :: ud(grid%ifirstxy:grid%ilastxy,plevp)         ! zonal deviation of zonal velocity
    real(r8) :: vd(grid%ifirstxy:grid%ilastxy,plevp)         ! zonal deviation of meridional velocity
    real(r8) :: wd(grid%ifirstxy:grid%ilastxy,plevp)         ! zonal deviation of vertical velocity
    real(r8) :: thd(grid%ifirstxy:grid%ilastxy,plevp)        ! zonal deviation of pot. temperature

    real(r8) :: vthp(grid%ifirstxy:grid%ilastxy,plevp)       ! zonal deviation of zonal velocity
    real(r8) :: wthp(grid%ifirstxy:grid%ilastxy,plevp)       ! zonal deviation of meridional velocity
    real(r8) :: uvp(grid%ifirstxy:grid%ilastxy,plevp)        ! zonal deviation of vertical velocity
    real(r8) :: uwp(grid%ifirstxy:grid%ilastxy,plevp)        ! zonal deviation of pot. temperature

    real(r8) :: rdiv(plevp)
    
    integer  :: ip_gm1g(plon,grid%jfirstxy:grid%jlastxy)     ! contains level index-1 where blocked points begin
    integer  :: zm_cnt(plevp)                                ! counter
    integer  :: i,j,k
    integer  :: nlons

    logical  :: has_zm(plevp,grid%jfirstxy:grid%jlastxy)                   ! .true. the (z,y) point is a valid zonal mean 
    integer  :: ip_gm1(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy) ! contains level index-1 where blocked points begin

    real(r8) :: vth(plevp,grid%jfirstxy:grid%jlastxy)                ! VTH flux
    real(r8) :: uv(plevp,grid%jfirstxy:grid%jlastxy)                 ! UV flux
    real(r8) :: wth(plevp,grid%jfirstxy:grid%jlastxy)                ! WTH flux
    real(r8) :: uw(plevp,grid%jfirstxy:grid%jlastxy)                 ! UW flux
    real(r8) :: u2d(plevp,grid%jfirstxy:grid%jlastxy)                ! zonally averaged U
    real(r8) :: v2d(plevp,grid%jfirstxy:grid%jlastxy)                ! zonally averaged V
    real(r8) :: th2d(plevp,grid%jfirstxy:grid%jlastxy)               ! zonally averaged TH
    real(r8) :: w2d(plevp,grid%jfirstxy:grid%jlastxy)                ! zonally averaged W
    real(r8) :: thig(grid%ifirstxy:grid%ilastxy,plevp,grid%jfirstxy:grid%jlastxy) ! interpolated pot. temperature

    real(r8) :: tmp2(grid%ifirstxy:grid%ilastxy)  
    real(r8) :: tmp3(grid%ifirstxy:grid%ilastxy,plevp)  

    integer :: beglat, endlat                          ! begin,end latitude indicies
    integer :: beglon, endlon                          ! begin,end longitude indicies

    beglon = grid%ifirstxy
    endlon = grid%ilastxy
    beglat = grid%jfirstxy
    endlat = grid%jlastxy

!omp parallel do private (i,j,k,pexf,psurf)
lat_loop1 : &
    do j = beglat, endlat
       do k = 1, plev
          do i = beglon, endlon
!-------------------------------------------------------------
! Calculate pressure and Exner function
!-------------------------------------------------------------
             pm(i,k,j) = 0.5_r8 * ( pe(i,k,j) + pe(i,k+1,j) )
             pexf      = (ps0 / pm(i,k,j))**cappa
!-------------------------------------------------------------
! Convert virtual temperature to temperature and calculate potential temperature
!-------------------------------------------------------------
             th(i,k,j) = pt(i,j,k) / (1._r8 + zvir*h2o(i,j,k)) 
             th(i,k,j) = th(i,k,j) * pexf
!-------------------------------------------------------------
! Calculate vertical velocity
!-------------------------------------------------------------
             w(i,k,j)  = - hscale * omga(i,k,j) / pm(i,k,j)
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
#ifdef CTEM_DIAGS
    write(iulog,*) '===================================================='
    do j = beglat,endlat
       write(iulog,'(''iam,myidxy_x,myidxy_y,j = '',4i4)') grid%iam,grid%myidxy_x,grid%myidxy_y,j
       write(iulog,'(20i3)') ip_gm1(:,j)
    end do
    if( grid%myidxy_x == 0 ) then
       do j = beglat,endlat
          write(iulog,*) '===================================================='
          write(iulog,'(''iam,myidxy_x,myidxy_y,j = '',4i4)') grid%iam,grid%myidxy_x,grid%myidxy_y,j
          write(iulog,'(20i3)') ip_gm1g(:,j)
       end do
       write(iulog,*) '===================================================='
#else
#ifdef SPMD    
    if( grid%myidxy_x == 0 ) then
#endif
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

#ifdef CTEM_DIAGS
    if( grid%myidxy_y == 12 ) then
       write(iulog,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       write(iulog,'(''iam,myidxy_x,myidxy_y,j = '',4i4)') grid%iam,grid%myidxy_x,grid%myidxy_y,beglat
       write(iulog,*) 'has_zm'
       write(iulog,'(20l2)') has_zm(:,beglat)
       write(iulog,*) 'ip_gm1g'
       write(iulog,'(20i4)') ip_gm1g(:,beglat)
       write(iulog,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
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
!      Zonal velocity
!-------------------------------------------------------------
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           u3(beglon,1,j), ui(beglon,k) )
!-------------------------------------------------------------
!      Meridional velocity
!-------------------------------------------------------------
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           v3(beglon,1,j), vi(beglon,k) )
!-------------------------------------------------------------
!      Vertical velocity
!-------------------------------------------------------------
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           w(beglon,1,j), wi(beglon,k) )
!-------------------------------------------------------------
!      Pot. Temperature
!-------------------------------------------------------------
          call vertinterp( nlons, nlons, plev, pm(beglon,1,j), pinterp, &
                           th(beglon,1,j), thi(beglon,k) )
       end do
#ifdef CTEM_DIAGS
       if( j == endlat ) then
       write(iulog,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       write(iulog,'(''iam,myidxy_x,myidxy_y,j = '',4i4)') grid%iam,grid%myidxy_x,grid%myidxy_y,j
       write(iulog,*) 'iref_p'
       write(iulog,'(5g15.7)') iref_p(:)
       write(iulog,'(''pm(endlon,:,'',i2,'')'')') j
       write(iulog,'(5g15.7)') pm(endlon,:,j)
       write(iulog,'(''u3(endlon,:,'',i2,'')'')') j
       write(iulog,'(5g15.7)') u3(endlon,:,j)
       write(iulog,*) 'ui(endlon,:)'
       write(iulog,'(5g15.7)') ui(endlon,:)
       write(iulog,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       end if
#endif

!-------------------------------------------------------------
! Calculate zonal averages
!-------------------------------------------------------------
       do k = ip_b+1, plevp
          if( has_zm(k,j) ) then
             where( ip_gm1(beglon:endlon,j) < k )
                ui(beglon:endlon,k)  = 0._r8
                vi(beglon:endlon,k)  = 0._r8
                wi(beglon:endlon,k)  = 0._r8
                thi(beglon:endlon,k) = 0._r8
             endwhere
          end if
       end do

       call par_xsum( grid, ui, plevp, um )
       call par_xsum( grid, vi, plevp, vm )
       call par_xsum( grid, wi, plevp, wm )
       call par_xsum( grid, thi, plevp, thm )
#ifdef CTEM_DIAGS
       if( j == endlat .and. grid%myidxy_y == 12 ) then
          write(iulog,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
          write(iulog,'(''iam,myidxy_x,myidxy_y,j = '',4i4)') grid%iam,grid%myidxy_x,grid%myidxy_y,j
          write(iulog,*) 'um after par_xsum'
          write(iulog,'(5g15.7)') um(:)
          write(iulog,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       end if
#endif
       do k = 1, ip_b
          um(k)     = um(k) * rplon
          vm(k)     = vm(k) * rplon
          wm(k)     = wm(k) * rplon
          thm(k)    = thm(k) * rplon
          u2d(k,j)  = um(k)
          v2d(k,j)  = vm(k)
          th2d(k,j) = thm(k)
          w2d(k,j)  = wm(k)
       end do
       do k = ip_b+1, plevp
          if( has_zm(k,j) ) then
             rdiv(k)   = 1._r8/count( ip_gm1g(:,j) >= k )
             um(k)     = um(k) * rdiv(k)
             vm(k)     = vm(k) * rdiv(k)
             wm(k)     = wm(k) * rdiv(k)
             thm(k)    = thm(k) * rdiv(k)
             u2d(k,j)  = um(k)
             v2d(k,j)  = vm(k)
             th2d(k,j) = thm(k)
             w2d(k,j)  = wm(k)
          else
             u2d(k,j)  = navp
             v2d(k,j)  = navp
             th2d(k,j) = navp
             w2d(k,j)  = navp
          end if
       end do

!-------------------------------------------------------------
! Calculate zonal deviations
!-------------------------------------------------------------
       do k = 1, ip_b
          ud(beglon:endlon,k)  = ui(beglon:endlon,k)  - um(k)
          vd(beglon:endlon,k)  = vi(beglon:endlon,k)  - vm(k)
          wd(beglon:endlon,k)  = wi(beglon:endlon,k)  - wm(k)
          thd(beglon:endlon,k) = thi(beglon:endlon,k) - thm(k)
       end do

       do k = ip_b+1, plevp
          if( has_zm(k,j) ) then
             where( ip_gm1g(beglon:endlon,j) >= k )
                ud(beglon:endlon,k)  = ui(beglon:endlon,k) - um(k)
                vd(beglon:endlon,k)  = vi(beglon:endlon,k) - vm(k)
                wd(beglon:endlon,k)  = wi(beglon:endlon,k) - wm(k)
                thd(beglon:endlon,k) = thi(beglon:endlon,k) - thm(k)
             elsewhere
                ud(beglon:endlon,k)  = 0._r8
                vd(beglon:endlon,k)  = 0._r8
                wd(beglon:endlon,k)  = 0._r8
                thd(beglon:endlon,k) = 0._r8
             endwhere
          end if
       end do

!-------------------------------------------------------------
! Calculate fluxes
!-------------------------------------------------------------
       do k = 1, ip_b
          vthp(:,k) = vd(:,k) * thd(:,k)
          wthp(:,k) = wd(:,k) * thd(:,k)
          uwp(:,k)  = wd(:,k) * ud(:,k)
          uvp(:,k)  = vd(:,k) * ud(:,k)
       end do

       do k = ip_b+1, plevp
          if( has_zm(k,j) ) then
             vthp(:,k) = vd(:,k) * thd(:,k)
             wthp(:,k) = wd(:,k) * thd(:,k)
             uwp(:,k)  = wd(:,k) * ud(:,k)
             uvp(:,k)  = vd(:,k) * ud(:,k)
          else
             vthp(:,k) = 0._r8
             wthp(:,k) = 0._r8
             uwp(:,k)  = 0._r8
             uvp(:,k)  = 0._r8
          end if
       end do

#ifdef CTEM_DIAGS
       if( j == endlat .and. grid%myidxy_y == 12 ) then
          write(iulog,*) '#################################################'
          write(iulog,*) 'DIAGNOSTICS before par_xsum'
          write(iulog,'(''iam,myidxy_x,myidxy_y,j = '',4i4)') grid%iam,grid%myidxy_x,grid%myidxy_y,j
          write(iulog,*) 'has_zm'
          write(iulog,*) has_zm(:,j)
          write(iulog,*) 'rdiv'
          write(iulog,'(5g15.7)') rdiv(:)
          write(iulog,*) 'wm'
          write(iulog,'(5g15.7)') wm(:)
          write(iulog,*) 'um'
          write(iulog,'(5g15.7)') um(:)
          write(iulog,*) 'uw'
          write(iulog,'(5g15.7)') uw(:)
          write(iulog,*) '#################################################'
       end if
#endif
       call par_xsum( grid, vthp, plevp, vth(1,j) )
       call par_xsum( grid, wthp, plevp, wth(1,j) )
       call par_xsum( grid, uvp, plevp, uv(1,j) )
       call par_xsum( grid, uwp, plevp, uw(1,j) )
#ifdef CTEM_DIAGS
       if( j == endlat .and. grid%myidxy_y == 12 ) then
          write(iulog,*) '#################################################'
          write(iulog,'(''iam,myidxy_x,myidxy_y,j = '',4i4)') grid%iam,grid%myidxy_x,grid%myidxy_y,j
          write(iulog,*) 'uw after par_xsum'
          write(iulog,'(5g15.7)') uw(:,j)
          write(iulog,*) '#################################################'
       end if
#endif
       do k = 1, ip_b
          vth(k,j) = vth(k,j) * rplon
          wth(k,j) = wth(k,j) * rplon
          uw(k,j)  = uw(k,j) * rplon
          uv(k,j)  = uv(k,j) * rplon
       end do
       do k = ip_b+1, plevp
          if( has_zm(k,j) ) then
             vth(k,j) = vth(k,j) * rdiv(k)
             wth(k,j) = wth(k,j) * rdiv(k)
             uw(k,j)  = uw(k,j) * rdiv(k)
             uv(k,j)  = uv(k,j) * rdiv(k)
          else
             vth(k,j) = navp
             wth(k,j) = navp
             uw(k,j)  = navp
             uv(k,j)  = navp
          end if
       end do

       thig(:,:,j) = thi(:,:)
    end do lat_loop3

!-------------------------------------------------------------
! Do the output
!-------------------------------------------------------------
    latloop: do j = beglat,endlat

!-------------------------------------------------------------
! zonal-mean output
!-------------------------------------------------------------
       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = vth(k,j)
       enddo
       call outfld( 'VTHzm', tmp3(grid%ifirstxy,:), 1, j )

       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = wth(k,j)
       enddo
       call outfld( 'WTHzm', tmp3(grid%ifirstxy,:), 1, j )

       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = uv(k,j)
       enddo
       call outfld( 'UVzm', tmp3(grid%ifirstxy,:), 1, j )
 
       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = uw(k,j)
       enddo
       call outfld( 'UWzm', tmp3(grid%ifirstxy,:), 1, j )
       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = u2d(k,j)
       enddo
       call outfld( 'Uzm', tmp3(grid%ifirstxy,:), 1, j )
       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = v2d(k,j)
       enddo
       call outfld( 'Vzm', tmp3(grid%ifirstxy,:), 1, j )
       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = w2d(k,j)
       enddo
       call outfld( 'Wzm', tmp3(grid%ifirstxy,:), 1, j )
       do k = 1,plevp
          tmp3(grid%ifirstxy,k) = th2d(k,j)
       enddo
       call outfld( 'THzm', tmp3(grid%ifirstxy,:), 1, j )

!-------------------------------------------------------------
! 3D output
!-------------------------------------------------------------
       do k = 1,plevp
          do i = beglon,endlon
             tmp3(i,k) = thig(i,k,j)
          enddo
       enddo
       call outfld( 'TH', tmp3, nlons, j )

!-------------------------------------------------------------
! horizontal output
!-------------------------------------------------------------
       tmp2(beglon:endlon) = ip_gm1(beglon:endlon,j)
       call outfld( 'MSKtem', tmp2, nlons, j  )

    enddo latloop

  end subroutine ctem_diags

!=================================================================================

  subroutine ctem_init()

  use hycoef,       only : hyai, hybi, ps0
  use phys_control, only : phys_getopts

!-------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------
    integer :: k
    logical :: history_waccm

    if (.not.do_circulation_diags) return

    rplon    = 1._r8/plon
    zm_limit = plon/3

!-------------------------------------------------------------
! Calculate reference pressure
!-------------------------------------------------------------
    do k = 1, plevp
       iref_p(k) = (hyai(k) + hybi(k)) * ps0
    end do
    if( masterproc ) then
       write(iulog,*) 'ctem_inti: iref_p'
       write(iulog,'(1p5g15.7)') iref_p(:)
    end if

!-------------------------------------------------------------
! Find level where hybrid levels become purely pressure 
!-------------------------------------------------------------
    ip_b = -1
    do k = 1,plev
       if( hybi(k) == 0._r8 ) ip_b = k
    end do

    call phys_getopts( history_waccm_out = history_waccm )

!-------------------------------------------------------------
! Initialize output buffer
!-------------------------------------------------------------
    call addfld ('VTHzm',(/ 'ilev' /),'A','MK/S','Meridional Heat Flux: 3D zon. mean', gridname='fv_centers_zonal' )
    call addfld ('WTHzm',(/ 'ilev' /),'A','MK/S','Vertical Heat Flux: 3D zon. mean', gridname='fv_centers_zonal' )
    call addfld ('UVzm', (/ 'ilev' /),'A','M2/S2','Meridional Flux of Zonal Momentum: 3D zon. mean', gridname='fv_centers_zonal' )
    call addfld ('UWzm', (/ 'ilev' /),'A','M2/S2','Vertical Flux of Zonal Momentum: 3D zon. mean', gridname='fv_centers_zonal' )

    call addfld ('Uzm',  (/ 'ilev' /),'A','M/S','Zonal-Mean zonal wind - defined on ilev', gridname='fv_centers_zonal' )
    call addfld ('Vzm',  (/ 'ilev' /),'A','M/S','Zonal-Mean meridional wind - defined on ilev', gridname='fv_centers_zonal' )
    call addfld ('Wzm',  (/ 'ilev' /),'A','M/S','Zonal-Mean vertical wind - defined on ilev', gridname='fv_centers_zonal' )
    call addfld ('THzm', (/ 'ilev' /),'A',  'K','Zonal-Mean potential temp - defined on ilev', gridname='fv_centers_zonal' )

    call addfld ('TH',   (/ 'ilev' /),'A','K',  'Potential Temperature', gridname='fv_centers' )
    call addfld ('MSKtem',horiz_only, 'A','1',  'TEM mask', gridname='fv_centers' )
    
!-------------------------------------------------------------
! primary tapes: 3D fields
!-------------------------------------------------------------
    call add_default ('VTHzm', 1, ' ')
    call add_default ('WTHzm', 1, ' ')
    call add_default ('UVzm' , 1, ' ')
    call add_default ('UWzm' , 1, ' ')
    call add_default ('TH' , 1, ' ')
    call add_default ('MSKtem',1, ' ')

    if (history_waccm) then
       call add_default ('MSKtem',7, ' ')
       call add_default ('VTHzm', 7, ' ')
       call add_default ('UVzm', 7, ' ')
       call add_default ('UWzm', 7, ' ')
       call add_default ('Uzm', 7, ' ')
       call add_default ('Vzm', 7, ' ')
       call add_default ('Wzm', 7, ' ')
       call add_default ('THzm', 7, ' ')
    end if

    if (masterproc) then
       write(iulog,*) 'ctem_inti: do_circulation_diags = ',do_circulation_diags
    endif

  end subroutine ctem_init

!================================================================================

subroutine ctem_readnl(nlfile)

   use namelist_utils,     only: find_group_name
   use units,              only: getunit, freeunit
   use spmd_utils,         only: mpicom, mstrid=>masterprocid, mpi_logical

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
   
   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'ctem_readnl'

   namelist /circ_diag_nl/ do_circulation_diags
   !-----------------------------------------------------------------------------
    
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'circ_diag_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, circ_diag_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(do_circulation_diags, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: do_circulation_diags")

end subroutine ctem_readnl

end module ctem
