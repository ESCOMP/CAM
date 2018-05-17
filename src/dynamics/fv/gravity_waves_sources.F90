module gravity_waves_sources

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plev, plevp
  use hycoef, only         : hypi

  implicit none
  save
  private

  public :: gws_src_fnct

  ! ghosting added by Francis Vitt -- 7 July 2008
  !
  ! moved from waccm to fv, changed source of psurf_ref
  !   -- S Santos -- 10 Aug 2011

  contains


!=================================================================

    subroutine gws_src_fnct(grid, u3,v3,pt, q3, pe, frontgf, frontga)

      use dynamics_vars, only: t_fvdycore_grid
      use physconst,     only: zvir, cappa, aearth => rearth

      implicit none

! Input/Output arguments
      type (t_fvdycore_grid), intent(in) :: grid    ! grid for XY decomp
      real(r8), intent(in) :: u3(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)          ! zonal velocity
      real(r8), intent(in) :: v3(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)          ! meridional velocity
      real(r8), intent(in) :: pt(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,plev)          ! virtual temperature
      real(r8), intent(in) :: q3(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,plev)          ! water constituent
      real(r8), intent(in) :: pe(grid%ifirstxy:grid%ilastxy,plevp,grid%jfirstxy:grid%jlastxy)         ! interface pressure

      real(r8), intent(out) :: frontgf(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)    ! Frontogenesis function
      real(r8), intent(out) :: frontga(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)    ! Frontogenesis angle

! Locals
      real(r8) :: psurf_ref                        ! surface reference pressure

      real(r8) :: ptemp(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)       ! temperature
      real(r8) :: pm(grid%ifirstxy:grid%ilastxy   ,plev,grid%jfirstxy:grid%jlastxy)       ! mid-point pressure
      real(r8) :: pexf                                ! Exner function

      real(r8) :: pty(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! temperature meridional gradient
      real(r8) :: ptx(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! temperature zonal gradient

      real(r8) :: uy(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! U-wind meridional gradient
      real(r8) :: ux(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! U-wind zonal gradient

      real(r8) :: vy(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! V-wind meridional gradient
      real(r8) :: vx(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)         ! V-wind zonal gradient

      real(r8) :: ptg(grid%ifirstxy-1:grid%ilastxy+1,plev,grid%jfirstxy-1:grid%jlastxy+1)  ! temperature ghosted
      real(r8) ::  ug(grid%ifirstxy-1:grid%ilastxy+1,plev,grid%jfirstxy-1:grid%jlastxy+1)  ! U-wind ghosted
      real(r8) ::  vg(grid%ifirstxy-1:grid%ilastxy+1,plev,grid%jfirstxy-1:grid%jlastxy+1)  ! V-wind ghosted

      real(r8) :: tglat                               ! tangent-latitude
      integer  :: i,j,k
      integer  :: im, ip
      integer  :: beglatxy, endlatxy, beglonxy, endlonxy
!-----------------------------------------------------------------------------------------

      beglatxy = grid%jfirstxy
      endlatxy = grid%jlastxy
      beglonxy = grid%ifirstxy
      endlonxy = grid%ilastxy

      pty(:,:,:) = 0._r8
      ptx(:,:,:) = 0._r8
      uy(:,:,:) = 0._r8
      ux(:,:,:) = 0._r8
      vy(:,:,:) = 0._r8
      vx(:,:,:) = 0._r8
      frontgf(:,:,:) = 0._r8
      frontga(:,:,:) = 0._r8

      psurf_ref = hypi(plev+1)

      !$omp parallel do private (i,j,k,pexf)
      do j = beglatxy, endlatxy
         do k = 1, plev
            do i = beglonxy, endlonxy

               ! Calculate pressure and Exner function
               pm(i,k,j) = 0.5_r8 * ( pe(i,k,j) + pe(i,k+1,j) )
               pexf      = (psurf_ref / pm(i,k,j))**cappa

               ! Convert virtual temperature to temperature and calculate potential temperature
               ptemp(i,k,j) = pt(i,j,k) / (1._r8 + zvir*q3(i,j,k)) * pexf

            end do
         end do
      end do

      call ghost_array(grid, ptemp, ptg)
      call ghost_array(grid, u3, ug)
      call ghost_array(grid, v3, vg)

      !$omp parallel do private (i,j,k)
      do k=1, plev
         do i=beglonxy, endlonxy
            do j=beglatxy, endlatxy

               ! Pot. Temperature
               pty(i,k,j) = ( ptg(i,k,j+1) - ptg(i,k,j-1) ) / (2._r8 * grid%dp)
               pty(i,k,j) = pty(i,k,j) / aearth

               ! U-wind
               uy(i,k,j) = ( ug(i,k,j+1) - ug(i,k,j-1) ) / (2._r8 * grid%dp)
               uy(i,k,j) = uy(i,k,j) / aearth

               ! V-wind
               vy(i,k,j) = ( vg(i,k,j+1) - vg(i,k,j-1) ) / (2._r8 * grid%dp)
               vy(i,k,j) = vy(i,k,j) / aearth

            end do
         end do
      end do

      !++rrg use 1.e-3 floor on cosine terms in the denominator of frontgf

      !$omp parallel do private (i,j,k,im,ip)
      do k=1, plev
         do j=beglatxy, endlatxy
            do i=beglonxy, endlonxy

               im = i-1
               ip = i+1

               ! Pot. Temperature
               ptx(i,k,j) = ( ptg(ip,k,j) - ptg(im,k,j) ) / (2._r8 * grid%dl)
               ptx(i,k,j) = ptx(i,k,j) / (aearth * (grid%cosp(j)+1.e-3_r8))

               ! U-wind
               ux(i,k,j) = ( ug(ip,k,j) - ug(im,k,j) ) / (2._r8 *grid%dl)
               ux(i,k,j) = ux(i,k,j) / (aearth * (grid%cosp(j)+1.e-3_r8))

               ! V-wind
               vx(i,k,j) = ( vg(ip,k,j) - vg(im,k,j) ) / (2._r8 *grid%dl)
               vx(i,k,j) = vx(i,k,j) / (aearth * (grid%cosp(j)+1.e-3_r8))

            end do
         end do
      end do

      !$omp parallel do private (i,j,k, tglat)
      do j=beglatxy, endlatxy

         tglat = grid%sinp(j) / (grid%cosp(j)+1.e-3_r8)

         do k=1, plev
            do i=beglonxy, endlonxy

               frontgf(i,k,j) =                                                                          &
                    - ptx(i,k,j)**2._r8 * (ux(i,k,j) - v3(i,k,j) * tglat / aearth)                          &
                    - pty(i,k,j)**2._r8 * vy(i,k,j)                                                         &
                    - ptx(i,k,j) * pty(i,k,j) * ( vx(i,k,j) + uy(i,k,j) + u3(i,k,j) * tglat / aearth )

            end do
         end do

      end do

      !--rrg use 1.e-3 floor on cosine terms in the denominator of frontgf

      !$omp parallel do private (i,j,k)
      do j=beglatxy, endlatxy
         do k=1, plev
            do i=beglonxy, endlonxy
               frontga(i,k,j) = atan2 ( pty(i,k,j) , ptx(i,k,j) + 1.e-10_r8 )
            end do
         end do
      end do

      return

    end subroutine gws_src_fnct

    subroutine ghost_array(grid, x, xg)

      ! subroutine added by Francis Vitt -- 7 July 2008

#if defined( SPMD )
      use mod_comm,      only: mp_send3d, mp_recv3d
#endif
      use dynamics_vars, only: T_FVDYCORE_GRID

      implicit none

      ! Input/Output arguments
      type (T_FVDYCORE_GRID), intent(in) :: grid    ! grid for XY decomp
      real(r8), intent(in)  :: x(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy)          ! zonal velocity
      real(r8), intent(out) :: xg(grid%ifirstxy-1:grid%ilastxy+1,plev,grid%jfirstxy-1:grid%jlastxy+1)          ! zonal velocity

      ! local variables
      real(r8) :: north(grid%ifirstxy:grid%ilastxy,plev)
      real(r8) :: south(grid%ifirstxy:grid%ilastxy,plev)
      real(r8) :: east(plev,grid%jfirstxy:grid%jlastxy)
      real(r8) :: west(plev,grid%jfirstxy:grid%jlastxy)
      integer  :: im, jm, km, ifirstxy, ilastxy, jfirstxy, jlastxy, iam, myidxy_y, nprxy_x
      integer  :: itot, dest, src, j, k

      im = grid%im
      jm = grid%jm
      km = grid%km

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      iam      = grid%iam
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x
      itot = ilastxy-ifirstxy+1

      xg(ifirstxy:ilastxy,:,jfirstxy:jlastxy) = x(ifirstxy:ilastxy,:,jfirstxy:jlastxy)

#if defined( SPMD )

  ! north
      call mp_send3d( grid%commxy, iam-nprxy_x, iam+nprxy_x, im, km, jm,      &
                      ifirstxy, ilastxy, 1, km, jfirstxy, jlastxy,           &
                      ifirstxy, ilastxy, 1, km, jfirstxy, jfirstxy, x )
      call mp_recv3d( grid%commxy, iam+nprxy_x, im, jm, km,                   &
                      ifirstxy, ilastxy, 1, km, jlastxy+1, jlastxy+1,        &
                      ifirstxy, ilastxy, 1, km, jlastxy+1, jlastxy+1, north )

  ! south
      call mp_send3d( grid%commxy, iam+nprxy_x, iam-nprxy_x, im, km, jm,      &
                      ifirstxy, ilastxy, 1, km, jfirstxy, jlastxy,           &
                      ifirstxy, ilastxy, 1, km, jlastxy,  jlastxy, x )
      call mp_recv3d( grid%commxy, iam-nprxy_x, im, jm, km,                   &
                      ifirstxy, ilastxy, 1, km, jfirstxy-1, jfirstxy-1,        &
                      ifirstxy, ilastxy, 1, km, jfirstxy-1, jfirstxy-1, south )

#endif

      if (itot .ne. im) then
#if defined( SPMD )

   ! east

         dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         call mp_send3d( grid%commxy, dest, src, im, km, jm,                  &
                         ifirstxy, ilastxy, 1, km, jfirstxy, jlastxy,        &
                         ifirstxy, ifirstxy, 1, km, jfirstxy, jlastxy, x )
         call mp_recv3d( grid%commxy, src, im, km, jm,                        &
                         ilastxy+1, ilastxy+1, 1, km, jfirstxy, jlastxy,     &
                         ilastxy+1, ilastxy+1, 1, km, jfirstxy, jlastxy, east )

   ! west

         dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         call mp_send3d( grid%commxy, dest, src, im, km, jm,                  &
                         ifirstxy, ilastxy, 1, km, jfirstxy, jlastxy,          &
                         ilastxy,  ilastxy, 1, km, jfirstxy, jlastxy, x )
         call mp_recv3d( grid%commxy, src, im, km, jm,                        &
                         ifirstxy-1, ifirstxy-1, 1, km, jfirstxy, jlastxy, &
                         ifirstxy-1, ifirstxy-1, 1, km, jfirstxy, jlastxy, west )
#endif

      else
!$omp parallel do private(j, k)
         do k = 1,km
            do j=jfirstxy,jlastxy
               east(k,j) = x(1, k,j)
               west(k,j) = x(im,k,j)
            enddo
         enddo
      endif

      if ( jfirstxy == 1 ) then
         xg(ifirstxy:ilastxy,:,jfirstxy-1) = xg(ifirstxy:ilastxy,:,jfirstxy)
      else
         xg(ifirstxy:ilastxy,:,jfirstxy-1) = south
      endif

      if ( jlastxy == jm ) then
         xg(ifirstxy:ilastxy,:,jlastxy+1) = xg(ifirstxy:ilastxy,:,jlastxy)
      else
         xg(ifirstxy:ilastxy,:,jlastxy+1) = north
      endif

      xg(ifirstxy-1,:,jfirstxy:jlastxy) = west
      xg( ilastxy+1,:,jfirstxy:jlastxy) = east

    end subroutine ghost_array
!=================================================================

  end module gravity_waves_sources
