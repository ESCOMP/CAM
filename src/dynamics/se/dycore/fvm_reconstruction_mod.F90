!==================================================================================
! The subroutine reconstruction is called from both a horizontal 
!    threaded region and a nested region for horizontal and 
!    vertical threading.  if horz_num_threads != horz_num_threads*vert_num_threads
!    then the timers calls will generate segfault...  So the simple solution is 
!    to deactivate them by default.
!
#define FVM_TIMERS .FALSE.
!==================================================================================
  !MODULE FVM_RECONSTRUCTION_MOD--------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                  !
  ! This module contains everything  to do (ONLY) a CUBIC (3rd order) reconstruction  !
  !                                                                           !
  ! IMPORTANT: the implementation is done for a ncfl > 1, which is not working !
  !            but it works for ncfl=1                                        !
  !
  ! This module has been recoded for multi-tracer efficiency (May, 2014)
  !
  !---------------------------------------------------------------------------!
module fvm_reconstruction_mod

  use shr_kind_mod,           only: r8=>shr_kind_r8
  use control_mod,            only: north, south, east, west, neast, nwest, seast, swest
  use cam_abortutils,         only: endrun
  use perf_mod,               only: t_startf, t_stopf


  implicit none
  private
!    integer, parameter, private:: nh = nhr+(nhe-1) ! = 2 (nhr=2; nhe=1)
    ! = 3 (nhr=2; nhe=2)
  public :: reconstruction, recons_val_cart, extend_panel_interpolate
!reconstruction_gradient,
contains
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE RECONSTRUCTION------------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
  ! DESCRIPTION: controls the cubic (3rd order) reconstructions:                      !
  !                                                                                   !
  ! CALLS: fillhalo_cubic, reconstruction_cubic                                       !
  ! INPUT: fcube    ...  tracer values incl. the halo zone                            !
  !        fvm    ...  structure incl. tracer values aso                            !                                   !
  ! OUTPUT:recons   ...  has the reconstruction coefficients (5) for the 3rd order    !
  !                      reconstruction: dx, dy, dx^2, dy^2, dxdy                     !
  !-----------------------------------------------------------------------------------!
  subroutine reconstruction(fcube,nlev_in,k_in,recons,irecons,llimiter,ntrac_in,&
       nc,nhe,nhr,nhc,nht,ns,nh,&
       jx_min,jx_max,jy_min,jy_max,&
       cubeboundary,halo_interp_weight,ibase,&
       spherecentroid,&
       recons_metrics,recons_metrics_integral,&
       rot_matrix,centroid_stretch,&
       vertex_recons_weights,vtx_cart,&
       irecons_actual_in)
    implicit none
    !
    ! dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc)
    !
    integer, intent(in) :: irecons
    integer, intent(in) :: nlev_in, k_in
    integer, intent(in) :: ntrac_in,nc,nhe,nhr,nhc,nht,ns,nh,cubeboundary
    real (kind=r8), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev_in,ntrac_in), intent(inout) :: fcube
    real (kind=r8), dimension(irecons,1-nhe:nc+nhe,1-nhe:nc+nhe,ntrac_in), intent(out)   :: recons
    integer, intent(in) :: jx_min(3), jx_max(3), jy_min(3), jy_max(3)
    integer              , intent(in):: ibase(1-nh:nc+nh,1:nhr,2)
    real (kind=r8), intent(in):: halo_interp_weight(1:ns,1-nh:nc+nh,1:nhr,2)
    real (kind=r8), intent(in):: spherecentroid(irecons-1,1-nhe:nc+nhe,1-nhe:nc+nhe)
    real (kind=r8), intent(in):: recons_metrics(3,1-nhe:nc+nhe,1-nhe:nc+nhe)
    real (kind=r8), intent(in):: recons_metrics_integral(3,1-nhe:nc+nhe,1-nhe:nc+nhe)
    integer              , intent(in):: rot_matrix(2,2,1-nhc:nc+nhc,1-nhc:nc+nhc)
    real (kind=r8), intent(in):: centroid_stretch(7,1-nhe:nc+nhe,1-nhe:nc+nhe)
    real (kind=r8), intent(in):: vertex_recons_weights(4,1:irecons-1,1-nhe:nc+nhe,1-nhe:nc+nhe)
    real (kind=r8), intent(in):: vtx_cart(4,2,1-nhc:nc+nhc,1-nhc:nc+nhc)

    logical,           intent(in) :: llimiter(ntrac_in)
    integer, optional, intent(in) :: irecons_actual_in

    integer                       :: irecons_actual

    real (kind=r8), dimension(1-nht:nc+nht,1-nht:nc+nht,3) :: f

    integer :: i,j,in,h,itr,k
    integer,               dimension(2,3)                              :: jx,jy

    if (present(irecons_actual_in)) then 
       irecons_actual = irecons_actual_in
    else
       irecons_actual = irecons
    end if

    jx(1,1)=jx_min(1); jx(2,1)=jx_max(1)-1
    jx(1,2)=jx_min(2); jx(2,2)=jx_max(2)-1
    jx(1,3)=jx_min(3); jx(2,3)=jx_max(3)-1

    jy(1,1)=jy_min(1); jy(2,1)=jy_max(1)-1
    jy(1,2)=jy_min(2); jy(2,2)=jy_max(2)-1
    jy(1,3)=jy_min(3); jy(2,3)=jy_max(3)-1

    !
    ! Initialize recons
    !    
    call zero_non_existent_ghost_cell(recons,irecons,cubeboundary,nc,nhe,ntrac_in)
    if (irecons_actual>1) then
       if(FVM_TIMERS) call t_startf('FVM:reconstruction:part#1')
       if (nhe>0) then
          do itr=1,ntrac_in
             call extend_panel_interpolate(nc,nhc,nhr,nht,ns,nh,&
                  fcube(:,:,k_in,itr),cubeboundary,halo_interp_weight,ibase,f(:,:,1),f(:,:,2:3))             
             call get_gradients(f(:,:,:),jx,jy,irecons,recons(:,:,:,itr),&
                  rot_matrix,centroid_stretch,nc,nht,nhe,nhc,irecons_actual)
          end do
       else
          do itr=1,ntrac_in
             call extend_panel_interpolate(nc,nhc,nhr,nht,ns,nh,&
                  fcube(:,:,k_in,itr),cubeboundary,halo_interp_weight,ibase,f(:,:,1))
             call get_gradients(f(:,:,:),jx,jy,irecons,recons(:,:,:,itr),&
                  rot_matrix,centroid_stretch,nc,nht,nhe,nhc,irecons_actual)
          end do
       end if
       if(FVM_TIMERS) call t_stopf('FVM:reconstruction:part#1')
       if(FVM_TIMERS) call t_startf('FVM:reconstruction:part#2')

       call slope_limiter(nhe,nc,nhc,fcube,jx,jy,k_in,nlev_in,ntrac_in,irecons,recons,&
            spherecentroid(:,1-nhe:nc+nhe,1-nhe:nc+nhe),&
            recons_metrics,vertex_recons_weights,vtx_cart,irecons_actual,llimiter,&
            cubeboundary)
       if(FVM_TIMERS) call t_stopf('FVM:reconstruction:part#2')
    end if
    if(FVM_TIMERS) call t_startf('FVM:reconstruction:part#3')
    select case (irecons_actual)
    case(1)
      do in=1,3
        do j=jy(1,in),jy(2,in)
          do i=jx(1,in),jx(2,in)
            recons(1,i,j,1:ntrac_in)         = fcube(i,j,k_in,1:ntrac_in)
            recons(2:irecons,i,j,1:ntrac_in) = 0.0_r8
          end do
        end do
      end do
    case(3)
!      do j=1-nhe,nc+nhe
!        do i=1-nhe,nc+nhe
      do in=1,3
        do j=jy(1,in),jy(2,in)
          do i=jx(1,in),jx(2,in)
             recons(1,i,j,1:ntrac_in)         = fcube(i,j,k_in,1:ntrac_in) &
                  - recons(2,i,j,1:ntrac_in)*spherecentroid(1,i,j) &
                  - recons(3,i,j,1:ntrac_in)*spherecentroid(2,i,j)
             recons(2,i,j,1:ntrac_in)         = recons(2,i,j,1:ntrac_in)
             recons(3,i,j,1:ntrac_in)         = recons(3,i,j,1:ntrac_in)
             recons(4:irecons,i,j,1:ntrac_in) = 0.0_r8
          end do
        end do
      end do
    case(6)
      do itr=1,ntrac_in
        do in=1,3
          do j=jy(1,in),jy(2,in)
            do i=jx(1,in),jx(2,in)
              recons(1,i,j,itr)  = fcube(i,j,k_in,itr) &
                   - recons(2,i,j,itr)*spherecentroid(1,i,j) &
                   - recons(3,i,j,itr)*spherecentroid(2,i,j) &
                   + recons(4,i,j,itr)*recons_metrics_integral(1,i,j) &
                   + recons(5,i,j,itr)*recons_metrics_integral(2,i,j) &
                   + recons(6,i,j,itr)*recons_metrics_integral(3,i,j)
              recons(2,i,j,itr) = recons(2,i,j,itr)                 &
                   - recons(4,i,j,itr)*2.0_r8*spherecentroid(1,i,j) &
                   - recons(6,i,j,itr)       *spherecentroid(2,i,j)
              recons(3,i,j,itr) = recons(3,i,j,itr)                &
                   - recons(5,i,j,itr)*2.0_r8*spherecentroid(2,i,j) &
                   - recons(6,i,j,itr)*spherecentroid(1,i,j)
              !
              ! recons(i,j,4:6) already set in get_gradients
              !
            end do
          end do
        end do
      end do
    case default
      write(*,*) "irecons out of range in get_ceof", irecons
    end select
    if(FVM_TIMERS) call t_stopf('FVM:reconstruction:part#3')

    !          recons(a,b,3) * (centroid(a,b,1)**2 - centroid(a,b,3)) + &
    !          recons(a,b,4) * (centroid(a,b,2)**2 - centroid(a,b,4)) + &
    !          recons(a,b,5) * (centroid(a,b,1) * centroid(a,b,2) - centroid(a,b,5)) + &


    !  call debug_halo(fvm,fcubenew,fpanel)
    !  call debug_halo_recons(fvm,recons,recons_trunk)
    !  call print_which_case(fvm)
    !
    !  call debug_halo_neighbor       (fvm,fotherface,fotherpanel)
    !  call debug_halo_neighbor_recons(fvm,recons,recons_trunk)
  end subroutine reconstruction
  !END SUBROUTINE RECONSTRUCTION--------------------------------------------CE-for FVM!

  subroutine get_gradients(f,jx,jy,irecons,gradient,rot_matrix,centroid_stretch,nc,nht,nhe,nhc,irecons_actual)
    implicit none
    integer,                                                         intent(in)   :: irecons,nc,nht,nhe,nhc,irecons_actual
    real (kind=r8), dimension(1-nht:nc+nht,1-nht:nc+nht,3),          intent(in)   :: f
    real (kind=r8), dimension(irecons,1-nhe:nc+nhe,1-nhe:nc+nhe),    intent(inout):: gradient
    integer,               dimension(2,3),                           intent(in)   :: jx,jy
    integer              , dimension(2,2,1-nhc:nc+nhc,1-nhc:nc+nhc), intent(in)   :: rot_matrix
    real (kind=r8), dimension(7,1-nhe:nc+nhe,1-nhe:nc+nhe),          intent(in)   :: centroid_stretch

    integer                     :: i,j,in
    real (kind=r8), dimension(2):: g
    real (kind=r8)              :: sign
    character(len=128)          :: errormsg 
    

    select case (irecons_actual)
    case(3)
      in=1
      do j=jy(1,in),jy(2,in)
        do i=jx(1,in),jx(2,in)
          !
          ! df/dx: 4-th-order finite difference: (-f(i+2)+8f(i+1)-8f(i-1)+f(i-2))/12dx
          !
          gradient(2,i,j) = -f(i+2,j  ,in)+8.0_r8*f(i+1,j  ,in)-8.0_r8*f(i-1,j  ,in)+f(i-2,j  ,in)
          gradient(3,i,j) = -f(i  ,j+2,in)+8.0_r8*f(i  ,j+1,in)-8.0_r8*f(i  ,j-1,in)+f(i  ,j-2,in)
        end do
      end do
      do in=2,3
        do j=jy(1,in),jy(2,in)
          do i=jx(1,in),jx(2,in)
            g(1) = -f(i+2,j  ,in)+8.0_r8*f(i+1,j  ,in)-8.0_r8*f(i-1,j  ,in)+f(i-2,j  ,in)
            g(2) = -f(i  ,j+2,in)+8.0_r8*f(i  ,j+1,in)-8.0_r8*f(i  ,j-1,in)+f(i  ,j-2,in)
            gradient(2:3,i,j) = MATMUL(rot_matrix(:,:,i,j),g(:))
          end do
        end do
      end do
      gradient(2,:,:) = centroid_stretch(1,:,:)*gradient(2,:,:)
      gradient(3,:,:) = centroid_stretch(2,:,:)*gradient(3,:,:)
    case (6)
      in=1
      do j=jy(1,in),jy(2,in)
        do i=jx(1,in),jx(2,in)
          !
          ! df/dx: 4-th-order finite difference: (-f(i+2)+8f(i+1)-8f(i-1)+f(i-2))/12dx
          !
          gradient(2,i,j) = -f(i+2,j  ,in)+ 8.0_r8*f(i+1,j  ,in)                 - 8.0_r8*f(i-1,j  ,in)+f(i-2,j  ,in)
          gradient(3,i,j) = -f(i  ,j+2,in)+ 8.0_r8*f(i  ,j+1,in)                 - 8.0_r8*f(i  ,j-1,in)+f(i  ,j-2,in)
          !
          ! d2f/dx2:
          !
          gradient(4,i,j) = -f(i+2,j  ,in)+16.0_r8*f(i+1,j  ,in)-30.0_r8*f(i,j,in)+16.0_r8*f(i-1,j  ,in)-f(i-2,j  ,in)
          gradient(5,i,j) = -f(i  ,j+2,in)+16.0_r8*f(i  ,j+1,in)-30.0_r8*f(i,j,in)+16.0_r8*f(i  ,j-1,in)-f(i  ,j-2,in)

          gradient(6,i,j) =  f(i+1,j+1,in)-       f(i+1,j-1,in)                 -       f(i-1,j+1,in)+f(i-1,j-1,in)
          !
          ! "stretching factors
          !
          gradient(2,i,j) = centroid_stretch(1,i,j)*gradient(2,i,j)
          gradient(3,i,j) = centroid_stretch(2,i,j)*gradient(3,i,j)
          
          gradient(4,i,j) = centroid_stretch(3,i,j)*gradient(4,i,j)+centroid_stretch(6,i,j)*gradient(2,i,j)
          gradient(5,i,j) = centroid_stretch(4,i,j)*gradient(5,i,j)+centroid_stretch(7,i,j)*gradient(3,i,j)
          
          gradient(6,i,j) = centroid_stretch(5,i,j)*gradient(6,i,j)
        end do
      end do
      do in=2,3
        if (SUM(rot_matrix(:,:,jx(1,in),jy(1,in)))==0) then
          sign=-1
        else
          sign=1
        end if
        do j=jy(1,in),jy(2,in)
          do i=jx(1,in),jx(2,in)
            g(1) = -f(i+2,j  ,in)+8.0_r8*f(i+1,j  ,in)-8.0_r8*f(i-1,j  ,in)+f(i-2,j  ,in)
            g(2) = -f(i  ,j+2,in)+8.0_r8*f(i  ,j+1,in)-8.0_r8*f(i  ,j-1,in)+f(i  ,j-2,in)
            gradient(2:3,i,j) = MATMUL(rot_matrix(:,:,i,j),g(:))

            g(1) = -f(i+2,j  ,in)+16.0_r8*f(i+1,j  ,in)-30.0_r8*f(i,j,in)+16.0_r8*f(i-1,j  ,in)-f(i-2,j  ,in)
            g(2) = -f(i  ,j+2,in)+16.0_r8*f(i  ,j+1,in)-30.0_r8*f(i,j,in)+16.0_r8*f(i  ,j-1,in)-f(i  ,j-2,in)
            gradient(4:5,i,j) = MATMUL(ABS(rot_matrix(:,:,i,j)),g(:))

            gradient(6,i,j) =  sign*(f(i+1,j+1,in)-       f(i+1,j-1,in)                 -       f(i-1,j+1,in)+f(i-1,j-1,in))
            !
            ! "stretching factors
            !
            gradient(2,i,j) = centroid_stretch(1,i,j)*gradient(2,i,j)
            gradient(3,i,j) = centroid_stretch(2,i,j)*gradient(3,i,j)
            
            gradient(4,i,j) = centroid_stretch(3,i,j)*gradient(4,i,j)+centroid_stretch(6,i,j)*gradient(2,i,j)
            gradient(5,i,j) = centroid_stretch(4,i,j)*gradient(5,i,j)+centroid_stretch(7,i,j)*gradient(3,i,j)
            
            gradient(6,i,j) = centroid_stretch(5,i,j)*gradient(6,i,j)
          end do
        end do
      end do
    case default
       write(errormsg, *) irecons
      call endrun('ERROR: irecons out of range in slope_limiter'//errormsg)
    end select
  end subroutine get_gradients


  subroutine slope_limiter(nhe,nc,nhc,fcube,jx,jy,k,nlev,ntrac,irecons,recons,spherecentroid,recons_metrics,&
       vertex_recons_weights,vtx_cart,irecons_actual,llimiter,cubeboundary)
    implicit none
    integer                                                           , intent(in) :: irecons_actual,k,nlev,ntrac
    integer                                                           , intent(in) :: irecons,nhe,nc,nhc
    real (kind=r8), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,ntrac)   , intent(inout) :: fcube
    real (kind=r8), dimension(irecons,1-nhe:nc+nhe,1-nhe:nc+nhe,ntrac), intent(inout):: recons
    integer,               dimension(2,3)                             , intent(in) :: jx,jy
    real (kind=r8), dimension(irecons-1,1-nhe:nc+nhe,1-nhe:nc+nhe)    , intent(in) :: spherecentroid
    real (kind=r8), dimension(3,1-nhe:nc+nhe,1-nhe:nc+nhe)            , intent(in) :: recons_metrics
    real (kind=r8), dimension(4,1:irecons-1,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in) :: vertex_recons_weights
    real (kind=r8), dimension(4,2,1-nhc:nc+nhc,1-nhc:nc+nhc)          , intent(in) :: vtx_cart
    logical,        dimension(ntrac)                                  , intent(in) :: llimiter
    integer                                                           , intent(in) :: cubeboundary

    real (kind=r8):: minval_patch,maxval_patch
    real (kind=r8):: phi, min_val, max_val,disc
    real (kind=r8):: v1,v2,v3,v4,vx1,vx2,vx3,vx4,vy1,vy2,vy3,vy4,r2,r3,r4,r5,r6,scx,scy,dx,dy,ex1,ex2,f0,val
    real (kind=r8):: m1,m2,m3

    real (kind=r8):: min_phi
    real (kind=r8):: extrema(2), xminmax(2),yminmax(2),extrema_value(13)

    real(kind=r8) :: invtmp  ! temporary to pre-compute inverses
    integer       :: itmp1,itmp2,i,j,in,vertex,n,itr,h

    real (kind=r8), dimension(-1:1) :: minval_array, maxval_array
    real (kind=r8), parameter :: threshold = 1.0E-40_r8
    character(len=128)        :: errormsg
    integer :: im1,jm1,ip1,jp1
           !
       ! fill in non-existent (in physical space) corner values to simplify
       ! logic in limiter code (min/max operation)
       !
    do itr=1,ntrac
       if (.not. llimiter(itr)) cycle
       if (cubeboundary>4) then
          select case(cubeboundary)
          case (nwest)
             do h=1,nhe+1
                fcube(0,nc+h  ,k,itr) = fcube(1-h,nc  ,k,itr)
                fcube(1-h,nc+1,k,itr) = fcube(1  ,nc+h,k,itr)
             end do
          case (swest)
             do h=1,nhe+1
                fcube(1-h,0,k,itr) = fcube(1,1-h,k,itr)
                fcube(0,1-h,k,itr) = fcube(1-h,1,k,itr)
             end do
          case (seast)
             do h=1,nhe+1
                fcube(nc+h,0  ,k,itr) = fcube(nc,1-h,k,itr)
                fcube(nc+1,1-h,k,itr) = fcube(nc+h,1,k,itr)
             end do
          case (neast)
             do h=1,nhe+1
                fcube(nc+h,nc+1,k,itr) = fcube(nc,nc+h,k,itr)
                fcube(nc+1,nc+h,k,itr) = fcube(nc+h,nc,k,itr)
             end do
          end select
       end if
    end do

    select case (irecons_actual)
       !
       ! PLM limiter
       !
    case(3)
       do in=1,3
          do j=jy(1,in),jy(2,in)
             do i=jx(1,in),jx(2,in)
                do itr = 1, ntrac
                   if (.not. llimiter(itr)) cycle
                   !rck combined min/max and unrolled inner loop
                   !minval_patch = MINVAL(fcube(i-1:i+1,j-1:j+1))
                   !maxval_patch = MAXVAL(fcube(i-1:i+1,j-1:j+1))
                   !DIR$ SIMD
                   do itmp2=-1,+1
                      itmp1 = j+itmp2
                      minval_array(itmp2) = min(fcube(i-1,itmp1,k,itr),fcube(i,itmp1,k,itr),fcube(i+1,itmp1,k,itr))
                      maxval_array(itmp2) = max(fcube(i-1,itmp1,k,itr),fcube(i,itmp1,k,itr),fcube(i+1,itmp1,k,itr))
                   enddo
                   minval_patch = min(minval_array(-1),minval_array(0),minval_array(1))
                   maxval_patch = max(maxval_array(-1),maxval_array(0),maxval_array(1))
                   min_phi=1.0_r8
                   !
                   ! coordinate bounds (could be pre-computed!)
                   !
                   xminmax(1) = min(vtx_cart(1,1,i,j),vtx_cart(2,1,i,j),vtx_cart(3,1,i,j),vtx_cart(4,1,i,j))
                   xminmax(2) = max(vtx_cart(1,1,i,j),vtx_cart(2,1,i,j),vtx_cart(3,1,i,j),vtx_cart(4,1,i,j))
                   yminmax(1) = min(vtx_cart(1,2,i,j),vtx_cart(2,2,i,j),vtx_cart(3,2,i,j),vtx_cart(4,2,i,j))
                   yminmax(2) = max(vtx_cart(1,2,i,j),vtx_cart(2,2,i,j),vtx_cart(3,2,i,j),vtx_cart(4,2,i,j))
                   !rck restructured loop
                   !DIR$ SIMD
                   do vertex=1,4
                      call recons_val_cart_plm(fcube(i,j,k,itr), vtx_cart(vertex,1,i,j), vtx_cart(vertex,2,i,j), spherecentroid(:,i,j), &
                           recons(1:3,i,j,itr), extrema_value(vertex))
                   end do
                   max_val = MAXVAL(extrema_value(1:4))
                   min_val = MINVAL(extrema_value(1:4))
                   if (max_val>maxval_patch.and.abs(max_val-fcube(i,j,k,itr))>threshold) then
                      phi = (maxval_patch-fcube(i,j,k,itr))/(max_val-fcube(i,j,k,itr))
                      if (phi<min_phi) min_phi=phi
                   end if
                   if (min_val<minval_patch.and.abs(min_val-fcube(i,j,k,itr))>threshold) then
                      phi = (minval_patch-fcube(i,j,k,itr))/(min_val-fcube(i,j,k,itr))
                      if (phi<min_phi) min_phi=phi
                   end if
                   ! Apply monotone limiter to all reconstruction coefficients
                   recons(2:3,i,j,itr)=min_phi*recons(2:3,i,j,itr)
                end do
             end do
          end do
       end do
       !
       ! PPM limiter (optimized)
       !
    case(6)
       !
       ! default branch
       !
       do in=1,3
          do j=jy(1,in),jy(2,in)
             do i=jx(1,in),jx(2,in)
                !
                ! coordinate bounds (could be pre-computed!)
                !
                vx1 = vtx_cart(1,1,i,j);  vy1 = vtx_cart(1,2,i,j)
                vx2 = vtx_cart(2,1,i,j);  vy2 = vtx_cart(2,2,i,j)
                vx3 = vtx_cart(3,1,i,j);  vy3 = vtx_cart(3,2,i,j)
                vx4 = vtx_cart(4,1,i,j);  vy4 = vtx_cart(4,2,i,j)
                xminmax(1) = min(vx1,vx2,vx3,vx4)
                xminmax(2) = max(vx1,vx2,vx3,vx4)
                yminmax(1) = min(vy1,vy2,vy3,vy4)
                yminmax(2) = max(vy1,vy2,vy3,vy4)
                scx = spherecentroid(1,i,j)
                scy = spherecentroid(2,i,j)
                m1 = recons_metrics(1,i,j)
                m2 = recons_metrics(2,i,j)
                m3 = recons_metrics(3,i,j)
                im1 = i-1; ip1 = i+1
                jm1 = j-1; jp1 = j+1
                do itr = 1, ntrac
                   if (.not. llimiter(itr)) cycle
                   !rck combined min/max and unrolled inner loop
                   !minval_patch = MINVAL(fcube(i-1:i+1,j-1:j+1))
                   !maxval_patch = MAXVAL(fcube(i-1:i+1,j-1:j+1))
                   minval_patch = min( min( fcube(im1,jm1,k,itr), fcube(i ,jm1,k,itr), fcube(ip1,jm1,k,itr) ), &
                                       min( fcube(im1,j  ,k,itr), fcube(i ,j  ,k,itr), fcube(ip1,j  ,k,itr) ), &
                                       min( fcube(im1,jp1,k,itr), fcube(i ,jp1,k,itr), fcube(ip1,jp1,k,itr) ) )
                   maxval_patch = max( max( fcube(im1,jm1,k,itr), fcube(i ,jm1,k,itr), fcube(ip1,jm1,k,itr) ), &
                                       max( fcube(im1,j  ,k,itr), fcube(i ,j  ,k,itr), fcube(ip1,j  ,k,itr) ), &
                                       max( fcube(im1,jp1,k,itr), fcube(i ,jp1,k,itr), fcube(ip1,jp1,k,itr) ) )
                   min_phi=1.0_r8

                   f0  = fcube(i,j,k,itr)
                   min_val = f0; max_val = f0!initialize min/max
                   !
                   ! compute min/max value at cell corners
                   !DIR$ SIMD
                   do vertex=1,4
                      val = f0
                      do itmp1=1,irecons-1
                         val = val + recons(itmp1+1,i,j,itr)*vertex_recons_weights(vertex,itmp1,i,j)
                      enddo
                      min_val = min(min_val,val)
                      max_val = max(max_val,val)
                   enddo
                   r2 = recons(2,i,j,itr);  r3 = recons(3,i,j,itr)
                   r4 = recons(4,i,j,itr);  r5 = recons(5,i,j,itr);  r6 = recons(6,i,j,itr)
                   ! Check if the quadratic is minimized within the element
                   ! Extrema in the interior of the element (there might be just one candidate)
                   ! DO NOT NEED ABS here, if disc<0 we have a saddle point (no maximum or minimum)
                   disc = 4.0_r8 * r4 *r5 - r6*r6
                   if (abs(disc) > threshold) then
                      ex1 = (r6*r3 - 2.0_r8*r5*r2)
                      ex2 = (r6*r2 - 2.0_r8*r4*r3)
                      disc = 1 /disc
                      ex1 = ex1*disc+scx
                      ex2 = ex2*disc+scy
                     if ( ex1 > xminmax(1)-threshold .and. ex1 < xminmax(2)+threshold .and. &
                           ex2 > yminmax(1)-threshold .and. ex2 < yminmax(2)+threshold ) then
                         dx = ex1 - scx; dy = ex2 - scy
                         v1 = f0 + r2*dx + r3*dy + r4*(m1+dx*dx) + r5*(m2+dy*dy) + r6*(m3+dx*dy)
                         max_val = max(max_val, v1)
                         min_val = min(min_val, v1)
                      end if
                   endif
                   !
                   ! Check all potential minimizer points along element boundaries
                   !

                   !
                   ! Top/bottom edge, y=const., du/dx=0
                   !
                   if (abs(r4) > threshold) then
                      invtmp = 1.0_r8 / (2.0_r8 * r4)
                      do n = 1,2
                         ex1 = scx+invtmp * (-r2 - r6 * (yminmax(n) - scy))
                         if ((ex1 > xminmax(1)-threshold) .and. (ex1 < xminmax(2)+threshold)) then
                            dx = ex1 - scx; dy = yminmax(n) - scy
                            v1 = f0 + r2*dx + r3*dy + r4*(m1+dx*dx) + r5*(m2+dy*dy) + r6*(m3+dx*dy)
                            max_val = max(max_val, v1)
                            min_val = min(min_val, v1)
                         endif
                      enddo
                   endif
                   !
                   ! Left/right edge, x=const., du/dy=0
                   !
                   if (abs(r5) > threshold) then
                      invtmp = 1.0_r8 / (2.0_r8 * r5)
                      do n = 1,2
                         ex1 = scy+invtmp * (-r3 - r6 * (xminmax(n) - scx))
                         if ((ex1 > yminmax(1)-threshold) .and. (ex1 < yminmax(2)+threshold)) then
                            dx = xminmax(n) - scx; dy = ex1 - scy
                            v1 = f0 + r2*dx + r3*dy + r4*(m1+dx*dx) + r5*(m2+dy*dy) + r6*(m3+dx*dy)
                            max_val = max(max_val, v1)
                            min_val = min(min_val, v1)
                         endif
                      enddo
                   endif
                   !
                   if (max_val>maxval_patch.and.abs(max_val-fcube(i,j,k,itr))>threshold) then
                      phi = (maxval_patch-fcube(i,j,k,itr))/(max_val-fcube(i,j,k,itr))
                      if (phi<min_phi) min_phi=phi
                   end if
                   if (min_val<minval_patch.and.abs(min_val-fcube(i,j,k,itr))>threshold) then
                      phi = (minval_patch-fcube(i,j,k,itr))/(min_val-fcube(i,j,k,itr))
                      if (phi<min_phi) min_phi=phi
                   end if
                   recons(2:6,i,j,itr)=min_phi*recons(2:6,i,j,itr)
                end do
             end do
          end do
       end do
    case default
       write(errormsg, *) irecons
      call endrun('ERROR: irecons out of range in slope_limiter'//errormsg)
    end select

  end subroutine slope_limiter

  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE RECONS_VAL_CART-----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
  ! DESCRIPTION: returns the value from the reconstruction (3rd order Taylor polynom) !
  !              at the point (cartx,carty) -> in cube CARTESIAN coordinates          !
  !                                                                                   !
  ! INPUT: fcube  ...  tracer values incl. the halo zone                              !
  !        cartx ...  x cartesian coordinate of the evaluation point                  !
  !        carty ...  y cartesian coordinate of the evaluation point                  !
  !        centroid..  x,y,x^2,y^2,xy                                                 !
  !        recons ...  array of reconstructed coefficients                            !
  ! OUTPUT: value ... evaluation at a given point                                     !
  !-----------------------------------------------------------------------------------!
  !DIR$ ATTRIBUTES FORCEINLINE :: recons_val_cart
  subroutine recons_val_cart(fcube, cartx, carty, centroid, pre_computed_metrics, recons, value)
    implicit none
    real(kind=r8), intent(in) :: fcube
    real(kind=r8), intent(in) :: cartx, carty
    real(kind=r8), dimension(1:5), intent(in) :: centroid
    real(kind=r8), dimension(3),   intent(in) :: pre_computed_metrics
    real(kind=r8), dimension(1:6), intent(in) :: recons
    real(kind=r8), intent(out) :: value
    real(kind=r8) :: dx, dy
    dx = cartx - centroid(1)
    dy = carty - centroid(2)
    ! Evaluate constant order terms
    value = fcube + &
         ! Evaluate linear order terms
         recons(2) * dx + &
         recons(3) * dy + &
         ! Evaluate second order terms
         recons(4) * (pre_computed_metrics(1) + dx*dx) + &
         recons(5) * (pre_computed_metrics(2) + dy*dy) + &
         recons(6) * (pre_computed_metrics(3) + dx*dy)
  END subroutine recons_val_cart
    !DIR$ ATTRIBUTES FORCEINLINE :: recons_val_cart_plm
    subroutine recons_val_cart_plm(fcube, cartx, carty, centroid, recons, value)
    implicit none
    real(kind=r8), intent(in) :: fcube
    real(kind=r8), intent(in) :: cartx, carty
    real(kind=r8), dimension(1:5), intent(in) :: centroid
    real(kind=r8), dimension(1:3), intent(in) :: recons
    real(kind=r8), intent(out) :: value
    real(kind=r8) :: dx, dy
    dx = cartx - centroid(1)
    dy = carty - centroid(2)
    ! Evaluate constant order terms
    value = fcube + &
         ! Evaluate linear order terms
         recons(2) * dx + &
         recons(3) * dy 
  END subroutine recons_val_cart_plm


  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE SLOPELIMITER_VAL----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
  ! DESCRIPTION: returns the value from the reconstruction (3rd order Taylor polynom) !
  !              at the point (cartx,carty) -> in cube CARTESIAN coordinates          !
  !                                                                                   !
  ! INPUT: value  ...  point value (calculated here by recons_val_cart)               !
  !        cell_value ...  tracer value (in the cell center) of the cell              !
  !        local_min ...  minmal value in the patch                                   !
  !        local_max ...  maximal value in the patch                                  !
  ! INPUT/OUTPUT: min_phi ... slope limiter, inout because we go through any possible !
  !                           extrema on the cell                                     !
  !-----------------------------------------------------------------------------------!
  subroutine slopelimiter_val(value, cell_value, local_min, local_max, min_phi)
    implicit none
    real (kind=r8), intent(in)    :: value, cell_value
    real (kind=r8), intent(in)    :: local_min, local_max
    real (kind=r8), intent(inout) :: min_phi
    real (kind=r8) :: phi

    ! Check against the minimum bound on the reconstruction
    if (value - cell_value > 1.0e-12_r8 * value) then
      phi = (local_max - cell_value) / (value - cell_value)
      if (phi < min_phi) then
        min_phi = phi
      endif
      ! Check against the maximum bound on the reconstruction
    elseif (value - cell_value < -1.0e-12_r8 * value) then
      phi = (local_min - cell_value) / (value - cell_value)
      if(phi < min_phi) then
        min_phi = phi
      endif
    endif
  end subroutine slopelimiter_val
  !END SUBROUTINE SLOPELIMITER_VAL------------------------------------------CE-for FVM!
  !DIR$ ATTRIBUTES FORCEINLINE :: dotproduct
  function dotproduct(w,f,ns)
    implicit none
    real (kind=r8)                          :: dotproduct
    real (kind=r8),dimension(:), intent(in) :: w,f      !dimension(ns)
    integer,                     intent(in) :: ns
    integer                                 :: k

    if(ns==3) then
      dotproduct = DotProduct_3(w,f)
    else
      dotproduct = DotProduct_gen(w,f,ns)
    endif

  end function dotproduct

  !DIR$ ATTRIBUTES FORCEINLINE :: DotProduct_gen
  function DotProduct_gen(w,f,ns)
    implicit none
    real (kind=r8)                          :: DotProduct_gen
    real (kind=r8),dimension(:), intent(in) :: w,f      !dimension(ns)
    integer,                     intent(in) :: ns
    integer                                 :: k
    DotProduct_gen = 0.0_r8
    do k=1,ns
       DotProduct_gen = DotProduct_gen+w(k)*f(k)
    end do
  end function DotProduct_gen

  ! special hard-coded version of the function where ns=3
  ! for performance optimization
  !DIR$ ATTRIBUTES FORCEINLINE :: DotProduct_3
  function DotProduct_3(w, f)
    IMPLICIT NONE
    REAL(KIND=r8), dimension(3), intent(in) :: w
    REAL(KIND=r8), dimension(3), intent(in) :: f
    REAL(KIND=r8) :: DotProduct_3
    DotProduct_3 = w(1)*f(1) + w(2)*f(2) + w(3)*f(3)
  end function DotProduct_3

  subroutine extend_panel_interpolate(nc,nhc,nhr,nht,ns,nh,fcube,cubeboundary,hWeight,ibase,&
       fpanel,fotherpanel)
    implicit none
    integer, intent(in) :: cubeboundary,nc,nhr,nht,nh,nhc,ns
    real (kind=r8),   &
         dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in)          :: fcube

    real (kind=r8), intent(in) :: hWeight(1:ns,1-nh:nc+nh,1:nhr,2)
    integer              , intent(in) :: ibase(1-nh:nc+nh,1:nhr,2)

    real (kind=r8)  , dimension(1-nht:nc+nht, 1-nht:nc+nht ), intent(out)           :: fpanel
    real   (kind=r8), dimension(1-nht:nc+nht,1-nht:nc+nht,2), intent(out), optional :: fotherpanel

    integer :: i, halo,ibaseref

    real (kind=r8), dimension(1-nhc:nc+nhc) :: ftmp
    !
    !  fpanel = 1.0E19 !dbg
    !
    !
    ! Stencil for reconstruction is:
    !
    !     ---------------------
    !     |   |   | i |   |   |
    !     ---------------------
    !     |   | i | i | i |   |
    !     ---------------------
    !     | i | i | R | i | i |
    !     ---------------------
    !     |   | i | i | i |   |
    !     ---------------------
    !     |   |   | i |   |   |
    !     ---------------------
    !
    ! where
    !
    !   "R" is cell for which we whish to do the reconstruction
    !   "i" is the stencil
    !
    !
    ! If one or more point in the stencil is on another panel(s) then we need to interpolate
    ! to a stencil that is an extension of the panel on which R is located
    ! (this is done using one dimensional cubic Lagrange interpolation along panel side)
    !
    ! Example: say that southern most "s" on Figure above is on another panels projection then the stencil becomes
    !
    !
    !   ---------------------------------
    !   |   |   |   |   |   | i |   |   |
    !   ----------------|----------------
    !   |   |   |   |   | i | i | i |   |
    !   ----------------|----------------
    !   |   |   |   | i | i | R | i | i |
    !   ----------------|----------------
    !   |   |   |   |   | i | i | i |   |
    !   ---------------------------------
    !   /   /   /   /   / S /S&i/ S / S /
    !  /---/---/---/---/---/---/---/---/
    ! /   /   /   /   /   /   /   /   /
    !/---/---/---/---/---/---/---/---/
    !
    !
    ! where "S" are the cell average values used for the cubic interpolation (located on the South panel)
    !
    !
    if (cubeboundary==0) then
      fpanel(1-nht:nc+nht,1-nht:nc+nht)=fcube(1-nht:nc+nht,1-nht:nc+nht)
    else if (cubeboundary==west) then
      !                                                       !
      !                                                       ! Case shown below: nhr=2, nhe=1, nht=nhr+nhe
      !                                                       ! (nhr = reconstruction width along x and y)
      !                                                       ! (nhe = max. Courant number)
      !                                                       !
      !                                                       !
      !    Figure below shows the element in question         ! In terms of data structure:
      !    (center element) and the surrounding elements      !
      !    on the element in question's projection            !     * "H" is on same panel average value
      !                                                       !     * "w" is west panel average values that need
      !    Notation: "0" marks the element boundaries         !       to be interpolated to main element
      !                                                       !       projection
      !    Elements to the west are on a different projection !     * "i" is extra halo required by the cubic
      !                                                       !       interpolation
      !     0                                                 !
      !     |0000                                             !
      !     |   |00000                                        !
      !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---x---x---x---x---x---x---x---x---x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   |   | i |   |   |   |   |   |   |   |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | i | H | H | H | H | H |   |   |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | w | H | H | H | H | H | H |   |   |
      !     |0000   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     |   |0000   0   |   |   |   0   |   |   |   0     !     |   | w | w | r | r | r | r | r | H | H |   |
      !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---00000000000000000---x---x---x---x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |\--|   \---0---------------0---------------0     !    -------------0---------------0---------------x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |\--|   \---0---------------0---------------0     !    -------------0---------------0---------------x
      !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |0000   |\--0---------------0---------------0     !    -------------0---------------0---------------x
      !     |   |0000   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---00000000000000000---x---x---x---x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w | r | r | r | r | r | H | H |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | w | H | H | H | H | H | H |   |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | i | H | H | H | H | H |   |   |   |
      !      0000   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !          0000   0   |   |   |   0   |   |   |   0     !     |   |   | i |   |   |   |   |   |   |   |   |
      !              000000000000000000000000000000000000     !    -x---x---x---x---x---x---x---x---x---x---x---x
      !
      !
      !      -2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8
      !
      !
      !
      ! fill in values (incl. halo) that are on the "main" panels projection
      !
      fpanel(1:nc+nht,1-nht:nc+nht)=fcube(1:nc+nht,1-nht:nc+nht)
      !
      ! fill in values that are on the west panels projection
      !
      do halo=1,nhr
        ftmp(:) = fcube(1-halo,:)   ! copy to a temporary
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref=ibase(i,halo,1)
          fpanel(1-halo ,i) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
        end do
      end do

      if (present(fotherpanel)) then
        !
        ! fill in values that are on the west panels projection
        !
        fotherpanel (1-nht:0,1-nht:nc+nht,1)=fcube(1-nht:0,1-nht:nc+nht)
        !
        do halo=1,nhr
          ftmp(:) = fcube(halo,:)   ! copy to a temporary
          do i=halo-nh,nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)
            !
            ! Exploit symmetry in interpolation weights
            !
            fotherpanel(halo,i,1)     = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
          end do
        end do
      end if
    else if (cubeboundary==east) then
      !
      ! north part is on different panel
      !
      ! stencil
      !
      ! CN<1 case                                             !
      !                                                       !
      !
      !
      !                                                 0     !
      !                                             0000|     !
      !                                         0000|   |     !
      !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---x---x---x---x---x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   |   |   |   |   |   |   | i |   |   |
      !     0---------------0---------------0--/    |--/|     !     x---------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   |   | H | H | H | H | H | i | i |   |
      !     0---------------0---------------0--/|   |--/|     !     x---------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   |   | H | H | H | H | H | H | e | i |   |
      !     0---------------0---------------0--/|   0000|     !     x---------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   0000|   |     !     |   | H | H | r | r | r | r | r | e | e |   |
      !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---00000000000000000---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     0---------------0---------------0--/|   |--/|     !     x---------------0---------------0---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     0---------------0---------------0--/|   |--/|     !     x---------------0---------------0---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     0---------------0---------------0--/|   0000|     !     x---------------0---------------0---x---x---x-
      !     0   |   |   |   0   |   |   |   0   0000|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---00000000000000000---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r | r | r | r | r | e | e |   |
      !     0---------------0---------------0--/|   |--/|     !     ----------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   | H | H | H | H | H | H | e | i |   |
      !     0---------------0---------------0--/|   |--/|     !     ----------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   |   |   | H | H | H | H | H | i | i |   |
      !     0---------------0---------------0--/|   0000      !     ----------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   0000          !     |   |   |   |   |   |   |   |   | i |   |   |
      !     000000000000000000000000000000000000              !     x---x---x---x---x---x---x---x---x---x---x---x-
      !
      !
      !      -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7
      !
      fpanel      (1-nht:nc     ,1-nht:nc+nht  )=fcube(1-nht:nc     ,1-nht:nc+nht)
      do halo=1,nhr
        ftmp(:) = fcube(nc+halo,:)   ! copy to a temporary
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref = ibase(i,halo,1)
          fpanel      (nc+halo   ,i  ) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
        end do
      end do

      if (present(fotherpanel)) then
        fotherpanel (nc+1 :nc+nht ,1-nht:nc+nht,1)=fcube(nc+1 :nc+nht ,1-nht:nc+nht) !
        do halo=1,nhr
          ftmp(:) = fcube(nc+1-halo,:)   ! copy to a temporary
          do i=halo-nh,nc+nh-(halo-1)
            ibaseref = ibase(i,halo,1)
            fotherpanel (nc+1-halo ,i,1) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
          end do
        end do
      end if

    else if (cubeboundary==north) then
      !
      ! north part is on different panel
      !
      ! stencil
      !
      ! CN<1 case
      !                                                      !   x---------------x---------------x---------------x
      !                                                      !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !0---\---\---\---0---\---\---\---0---\---\---\---0     !   x---------------x---------------x---------------x
      ! 0   \   \   \   0   \   \   \   0   \   \   \   0    !   |   | i | i | n | n | n | n | n | n | i | i |   |
      !  0---\---\---\---0---\---\---\---0---\---\---\---0   !   x---------------x---------------x---------------x
      !   0   \   \   \   0   \   \   \   0   \   \   \   0  !   | i | i | n | n | n | n | n | n | n | n | i | i |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r | r | r | r | r | r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   | H | H | H | H | H | H | H | H |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   | H | H | H | H | H | H |   |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---x---x---x---x---x---x---x---x---x
      !
      !    -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8    !    -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8
      !
      ! fill in values that are on the same projection as "main" element
      fpanel      (1-nht:nc+nht ,1-nht:nc)=fcube(1-nht:nc+nht ,1-nht:nc)
      do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref = ibase(i,halo,1)
          fpanel      (i,nc+halo    ) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,nc+halo),ns) !north
        end do
      end do
      if (present(fotherpanel)) then
        ! fill in halo for north element
        fotherpanel (1-nht:nc+nht ,nc+1:nc+nht,1)=fcube(1-nht:nc+nht ,nc+1:nc+nht)
        !
        do halo=1,nhr
          do i=halo-nh,nc+nh-(halo-1)
            ibaseref = ibase(i,halo,1)
            fotherpanel (i,nc+1-halo,1) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
          end do
        end do
      end if


    else if (cubeboundary==south) then
      !
      ! south part is on different panel
      !
      ! stencil
      !
      !                                                      !
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---x---x---x---x---x---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   | H | H | H | H | H | H |   |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   | H | H | H | H | H | H | H | H |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r | r | r | r | r | r | H | H |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   /   /   /   0   /   /   /   0   /   /   /   0  !   | i | i | s | s | s | s | s | s | s | s | i | i |
      !  0---/---/---/---0---/---/---/---0---/---/---/---0   !   x---------------x---------------x---------------x
      ! 0   /   /   /   0   /   /   /   0   /   /   /   0    !   |   | i | i | s | s | s | s | s | s | i | i |   |
      !0---/---/---/---0---/---/---/---0---/---/---/---0     !   x---------------x---------------x---------------x
      !                                                      !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !
      !     0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9           !     0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
      !
      ! fill in values that are on the same projection as "main" element (marked with "i" in Figure above)
      !
      fpanel      (1-nht:nc+nht,1:nc+nht  )=fcube(1-nht:nc+nht,1:nc+nht)
      do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref=ibase(i,halo,1)!ibase(i,halo,2)
          fpanel      (i,1-halo ) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,1-halo),ns)  !south
        end do
      end do
      if (present(fotherpanel)) then
        fotherpanel (1-nht:nc+nht,1-nht:0 ,1)=fcube(1-nht:nc+nht,1-nht:0 )
        do halo=1,nhr
          do i=halo-nh,nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)!ibase(i,halo,2)
            fotherpanel (i,  halo,1) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,  halo),ns)
          end do
        end do
      end if
    else if (cubeboundary==swest) then
      !
      ! south and west neighboring cells are on different panel
      !
      ! stencil
      !
      !
      ! CN<1 case
      !
      !
      !
      !     |000000000000000000000000000000000000   !   x---x---x---x---x---x---x---x---x---x---x---x---x
      !  0000   0   |   |   |   0   |   |   |   0   !   |   |   |   |   |   |   |   |   |   |   |   |   |
      ! 0   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   |   | w | H | H | H | H | H |   |   |   |
      ! |   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w | H | H | H | H | H | H |   |   |
      ! |   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w | r | r | r | r | r | H | H |   |
      ! |   |000000000000000000000000000000000000   !   x---x---x---x---00000000000000000---x---x---x---x
      ! |0000   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! 0   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
      ! |   |   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! | -/|   000000000000000000000000000000000   !   x---x---x---x---00000000000000000---x---x---x---x
      ! |/  | 0    /   /   /   0   /   /   /   0    !   |   |   | w |   | s | s | s | s | s | s |   |   |
      ! |   0-----/---/---/---0---/---/---/---0     !   x---------------x---------------x---------------x
      ! | 0      /   /   /   0   /   /   /   0      !   |   |   |   | s | s | s | s | s | s |   |   |   |
      ! 0-------/---/---/---0---/---/---/---0       !   x---------------x---------------x---------------x
      !
      !
      !  -1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |       !   |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
      !
      ! fill in values that are on the same projection as "main" element (marked with "i" in Figure above)
      !
      fpanel(1:nc+nht,1:nc+nht)=fcube(1:nc+nht,1:nc+nht)
      !
      ! fill in west part (marked with "w" on Figure above) and south part (marked with "s")
      !
      do halo=1,nhr
        ftmp(:)  = fcube(1-halo,:)   ! copy to a temporary
        do i=max(halo-nh,0),nc+nh-(halo-1)
          ibaseref=ibase(i,halo,1)!ibase(i,halo,1)
          fpanel(1-halo ,i) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)   !west
          fpanel(i,1-halo ) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,1-halo) ,ns)  !south
        end do
      end do
      !
      ! corner value
      !
      fpanel(0,0)  =0.25_r8*(fpanel(0,1)+fpanel(1,0)+fpanel(-1,0)+fpanel(0,-1))
      !
      ! ****************************************************************
      !
      ! fill halo for reconstruction on south neighbor panel projection
      !
      ! ****************************************************************
      !
      ! On the south panel projection the neighbors are arragened as follow (nwest case):
      !
      !
      ! \
      !  \    p
      !   \
      !    \-----
      !    |
      !  w |  s
      !    |
      !
      !
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   | p 0 p | p | p | p 0 p |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   | w | wp0 p | p | p | p 0 p | p |   |   |
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   |   | w | w | r | r | r | r | r | i | i |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   | w | i | i | i | i | i | i |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   | i | i | i | i | i |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   |   |   |   |   |   |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---
      !
      !
      ! fill values on same panel projection ("r" and "i" on Figure above)
      !
      if (present(fotherpanel)) then
        fotherpanel(1:nc+nht,1-nht:0,1)  = fcube(1:nc+nht,1-nht:0)
        !
        ! compute interpolated cell average values in "p" cells on Figure on above
        !
        do halo=1,nhr
          do i=max(halo-nh,0),nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)
            !
            ! use same weights as interpolation south from main panel (symmetric)
            !
            fotherpanel(i,halo,1)  = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,halo),ns)
          end do
        end do
        !
        ! compute interpolated cell average values in "w" cells on Figure on above
        !
        do halo=1,nhr
          do i=nc+halo-nhr,nc+1
            ibaseref=ibase(i,halo,2)-nc
            !
            ! fotherpanel indexing follows main panel indexing
            ! fcube indexing most be "rotated":
            !
            ! ===============================
            ! |              |              |
            ! |  W      ^    |   S          |
            ! |         |    |              |
            ! |       x |    |              |
            ! |         |    |              |
            ! !              |              |
            ! !   <-----     |              |
            ! !      y       |              |
            ! !              |              |
            ! ===============================
            !
            fotherpanel(1-halo,i-nc,1)  = dotproduct(hWeight(:,i,halo,2),fcube(ibaseref:ibaseref+ns-1,halo),ns)
          end do
        end do
        fotherpanel(0,1,1) = 0.25_r8*(fotherpanel(-1,1,1)+fotherpanel(1,1,1)+fotherpanel(0,2,1)+fotherpanel(0,0,1))
        !
        ! ****************************************************************
        !
        ! fill halo for reconstruction on west neighbor panel projection
        !
        ! ****************************************************************
        !
        ! On the west panel projection the neighbors are arragened as follow (seast case):
        !
        !   --------
        !   |      |
        !   |  w   |    p
        !   |      |
        !   -------\
        !           \
        !       s    \
        !
        !
        !
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   | i |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   | i | i | e |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   | i | i | r | e | e |   |   |   |   |   |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   |   | s | s | se| e |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   | s | s |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---
        !
        !
        ! fill values on same panel projection ("r" and "i" on Figure above)
        !
        fotherpanel(1-nht:nc,1:nc+nht,2)  = fcube(1-nht:nc,1:nc+nht)
        !
        ! compute interpolated cell average values in "p" cells on Figure on above
        !
        do halo=1,nhr
          ftmp(:) = fcube(halo,:)   ! copy to a temporary
          do i=max(halo-nh,0),nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)
            !
            ! use same weights as interpolation south from main panel (symmetric)
            !
            fotherpanel(halo,i,2)  = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        !
        ! compute interpolated cell average values in "s" cells on Figure on above
        !
        do halo=1,nhr
          ftmp(:) = fcube(halo,:)   ! copy to a temporary
          do i=nc+halo-nhr,nc+1
            ibaseref=ibase(i,halo,2)-nc
            !
            ! fotherpanel indexing follows main panel indexing
            ! fcube indexing most be "rotated":
            !
            ! ===============================
            ! |              |              |
            ! |  W      ^    |   S          |
            ! |         |    |              |
            ! |       x |    |              |
            ! |         |    |              |
            ! !              |              |
            ! !   <-----     |              |
            ! !      y       |              |
            ! !              |              |
            ! ===============================
            !
            fotherpanel(i-nc,1-halo,2)  = dotproduct(hWeight(:,i,halo,2),ftmp(ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        fotherpanel(1,0,2) = 0.25_r8*(fotherpanel(0,0,2)+fotherpanel(2,0,2)+fotherpanel(1,-1,2)+fotherpanel(1,1,2))
      end if
    else if (cubeboundary==seast) then
      !
      ! south and east neighboring cells are on different panel
      !
      !
      !
      ! 000000000000000000000000000000000000|
      ! 0   |   |   |   0   |   |   |   0   0000    !   |   |   |   |   |   |   |   |   |   |   |   |   |
      ! 0---------------0---------------0--\|   0   !   x---------------x---------------x---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   |   |   |   | H | H | H | H |   |   |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------x---------------x---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   |   |   | H | H | H | H | H | e |   |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------x---------------x---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r | r | r | r | r | e | e |   |   |
      ! 000000000000000000000000000000000000|   |   !   x---x---x---x---00000000000000000---x---x---x---x
      ! 0   |   |   |   0   |   |   |   0   0000|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 0---------------0---------------0--\|   0   !   x---------------0---------------0---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------0---------------0---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------0---------------0---------------x
      ! 0   |   |   |   0   |   |   |   0   |   |   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 000000000000000000000000000000000   |\- |   !   x---x---x---x---00000000000000000---x---x---x---x
      !  0   \   \   \   0   \   \   \    0 |  \|   !   |   |   | s | s | s | s | s | s |s/e| e |   |   |
      !   0---\---\---\---0---\---\---\-----0   |   !   x---------------x---------------x---------------x
      !    0   \   \   \   0   \   \   \      0 |   !   |   |   |   | s | s | s | s | s | s |   |   |   |
      !     0---\---\---\---0---\---\---\-------0   !   x---------------x---------------x---------------x
      !
      !
      fpanel       (1-nht:nc,1:nc+nht)=fcube(1-nht:nc,1:nc+nht)
      !
      ! east
      !
      do halo=1,nhr
        ftmp(:) = fcube(nc+halo,:)   ! copy to a temporary
        do i=max(halo-nh,0),nc+nh-(halo-1)
          ibaseref = ibase(i,halo,1)
          fpanel(nc+halo,i) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
        end do
      end do
      !
      ! south
      !
      do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
          ibaseref = ibase(i,halo,2)
          fpanel(i,1-halo ) = dotproduct(hWeight(:,i,halo,2),fcube(ibaseref:ibaseref+ns-1,1-halo),ns)  !south
        end do
      end do
      fpanel(nc+1,0   )=0.25_r8*(&
           fpanel(nc+1,1)+fpanel(nc,0)+fpanel(nc+2,0)+fpanel(nc+1,-1))
      !
      ! ****************************************************************
      !
      ! fill halo for reconstruction on south neighbor panel projection
      !
      ! ****************************************************************
      !
      ! On the south panel projection the neighbors are arragened as follow (neast case):
      !
      !
      !             /
      !       P    /
      !           /
      !    ------/
      !    |     |  E
      !    |  S  |
      !    |     |
      !
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   |   |   |   |   |   |   |   |   |
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   | n 0 n | n | n | n 0 n |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   | n | n 0 n | n | n | n 0 ne| e |   |   |
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   | i | i | r | r | r | r | r | e | e |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   | i | i | i | i | i | i | e |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   | i | i | i | i | i |   |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   |   |   |   |   |   |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      !
      !
      !
      ! fill values on same panel projection ("r" and "i" on Figure above)
      !
      if (present(fotherpanel)) then
        fotherpanel(1-nht:nc,1-nht:0,1)  = fcube(1-nht:nc,1-nht:0)
        !
        !
        ! fill in "n" on Figure above
        !
        do halo=1,nhr
          do i=halo-nh,min(nc+nh-(halo-1),nc+1)
            ibaseref = ibase(i,halo,2)
            fotherpanel (i,halo,1) = dotproduct(hWeight(:,i,halo,2),fcube(ibaseref:ibaseref+ns-1,  halo),ns)
          end do
        end do
        !
        ! fill in "e" on Figure above
        !
        do halo=1,nhr
          do i=0,nht-halo!nc+nh-(halo-1)
            ibaseref = ibase(i,halo,1)
            !
            ! fother panel follows indexing on main panel
            !
            ! use symmetry for weights (same weights as East from main panel but for south panel
            ! projection the indecies are rotated)
            !
            fotherpanel (nc+halo ,1-i,1) = dotproduct(hWeight(:,i,halo,1),fcube(nc+ibaseref:nc+ibaseref+ns-1,halo),ns)
          end do
        end do
        fotherpanel(nc+1,1,1) = 0.25_r8*(fotherpanel(nc+2,1,1)+fotherpanel(nc,1,1)&
             +fotherpanel(nc+1,2,1)+fotherpanel(nc+1,0,1))

        !
        ! ****************************************************************
        !
        ! fill halo for reconstruction on east neighbor panel projection
        !
        ! ****************************************************************
        !
        ! On the south panel projection the neighbors are arragened as follow (neast case):
        !
        !
        !             |     |
        !         P   |  E  |
        !             |-----|
        !            /
        !           /    S
        !          /
        !
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   | i |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   | w | i | i |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   | w | w | r | i | i |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   |   |   |   |   |   | w | ws| s | s |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   | s | s |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        !
        !
        !
        ! fill values on same panel projection ("r" and "i" on Figure above)
        !
        fotherpanel(nc+1:nc+nht,1:nc+nht,2)  = fcube(nc+1:nc+nht,1:nc+nht)
        !
        !
        ! fill in "w" on Figure above
        !
        do halo=1,nhr
          ftmp(:) = fcube(nc+1-halo,:)   ! copy to a temporary
          do i=0,nc+nh-(halo-1)
            ibaseref = ibase(i,halo,1)
            fotherpanel(nc+1-halo,i,2) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        !
        ! fill in "s" on Figure above
        !
        do halo=1,nhr
          ftmp(:) = fcube(nc+1-halo,:)   ! copy to a temporary
          do i=nc+1-nht+halo,nc+1
            !
            !
            ! !  P           |  E
            ! !              |
            ! !              |
            ! ================
            ! |              |
            ! |  S      |    |  <----- y
            ! |         |    |           ^
            ! |       x |    |           |
            ! |         v    |           |
            ! !              |           |
            ! !    ----->    |           x
            ! !      y       |
            ! !              |
            ! ================
            !
            !
            ! shift (since we are using south weights from main panel interpolation
            !
            ibaseref = ibase(i,halo,2)-nc
            !
            ! fotherpanel index: reverse
            !
            ! fcube index: due to rotation (see Figure above)
            !
            fotherpanel(nc+(nc+1-i),1-halo,2) = dotproduct(hWeight(:,i,halo,2),ftmp(ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        fotherpanel(nc,0,2) = 0.25_r8*(fotherpanel(nc+1,0,2)+fotherpanel(nc-1,0,2)&
             +fotherpanel(nc,1,2)+fotherpanel(nc,-1,2))
      end if
    else if (cubeboundary==nwest) then
      !
      !
      ! 0-------\---\---\---0---\---\---\---0       !   --------x---------------x---------------x
      ! | 0      \   \   \   0   \   \   \   0      !   |   | n | n | n | n | n | n |   |   |   |
      ! |   0-----\---\---\---0---\---\---\---0     !   --------x---------------x---------------x
      ! |   | 0    \   \   \   0   \   \   \   0    !   | w | a | n | n | n | n | n | n |   |   |
      ! |\  |   000000000000000000000000000000000   !   --------00000000000000000---------------x
      ! | -\|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |\--0---------------0---------------0   !   --------0---------------0---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |\--0---------------0---------------0   !   --------0---------------0---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! 0   |\--0---------------0---------------0   !   --------0---------------0---------------x
      ! |0000   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |000000000000000000000000000000000000   !   --------00000000000000000---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w | r | r | r | r | r | H | H |   |
      ! |   |\--0---------------0---------------0   !   --------x---------------x---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   |   | w | H | H | H | H | H | H |   |   |
      ! |   |\--0---------------0---------------0   !   --------x---------------x---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   |   |   | H | H | H | H | H |   |   |   |
      ! 0   |\--0---------------0---------------0   !   --------x---------------x---------------x
      !  0000   0   |   |   |   0   |   |   |   0   !   |   |   |   |   |   |   |   |   |   |   |
      !      000000000000000000000000000000000000   !   --------x---------------x---------------x
      !
      !
      !
      fpanel(1:nc+nht,1-nht:nc)=fcube(1:nc+nht,1-nht:nc)
      !
      ! west
      !
      do halo=1,nhr
        ftmp(:) = fcube(1-halo,:)   ! copy to a temporary
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
          ibaseref=ibase(i,halo,1)
          fpanel(1-halo ,i) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
        end do
      end do
      !
      ! north
      !
      do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref = ibase(i,halo,2)
           fpanel(i,nc+halo) = dotproduct(hWeight(:,i,halo,2),fcube(ibaseref:ibaseref+ns-1,nc+halo  ),ns) !north
         end do
       end do
       fpanel(0   ,nc+1)=0.25_r8*(&
            fpanel(0,nc)+fpanel(1,nc+1)+fpanel(-1,nc+1)+fpanel(0,nc+2))
       !
       ! ****************************************************************
       !
       ! fill halo for reconstruction on north neighbor panel projection
       !
       ! ****************************************************************
       !
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   |   |   |   |   |   |   |   |   |   |   |
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   |   |   | i | i | i | i | i |   |   |   |
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   |   | w | i | i | i | i | i | i |   |   |
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   | w | w | r | r | r | r | r | i | i |   |
       !x---x---x---x---00000000000000000---x---x---x---x
       !|   |   | w | ws0 s | s | s | s 0 s | s |   |   |
       !x---x---x---x---0---x---x---x---0---x---x---x---x
       !|   |   |   | s 0 s | s | s | s 0 s |   |   |   |
       !x---x---x---x---0---x---x---x---0---x---x---x---x
       !|   |   |   |   0   |   |   |   0   |   |   |   |
       !x---x---x---x---0---x---x---x---0---x---x---x---x
       !|   |   |   |   0   |   |   |   0   |   |   |   |
       !x---x---x---x---00000000000000000---x---x---x---x
       !
       !
       ! fill values on same panel projection ("r" and "i" on Figure above)
       !
       if (present(fotherpanel)) then
         fotherpanel(1:nc+nht,nc+1:nc+nht,1)  = fcube(1:nc+nht,nc+1:nc+nht)
         !
         !
         ! fill in "s" on Figure above
         !
         ! (use code from north above)
         !
         do halo=1,nhr
           do i=max(halo-nh,0),nc+nh-(halo-1)
             ibaseref = ibase(i,halo,2)
             fotherpanel(i,nc+1-halo,1) = dotproduct(hWeight(:,i,halo,2),fcube(ibaseref:ibaseref+ns-1,nc+1-halo  ),ns)
           end do
         end do
         !
         ! fill in "w" on Figure above
         !
         ! (use code from west above)
         !
         do halo=1,nhr
           do i=nc+1-nht+halo,nc+1
             ibaseref=ibase(i,halo,1)-nc
             fotherpanel(1-halo,nc-(i-(nc+1)),1) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
           end do
         end do
         fotherpanel(0,nc,1)=0.25_r8*(&
              fotherpanel(1,nc,1)+fotherpanel(-1,nc,1)+fotherpanel(0,nc+1,1)+fotherpanel(0,nc-1,1))

         !
         ! ****************************************************************
         !
         ! fill halo for reconstruction on west neighbor panel projection
         !
         ! ****************************************************************
         !
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   | n | n |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   | n | n | ne| e |   |   |   |   |   |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   | i | i | r | e | e |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   | i | i | e |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   | i |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---
         !
         !
         ! fill values on same panel projection ("r" and "i" on Figure above)
         !
         fotherpanel(1-nht:nc,1-nht:nc,2)  = fcube(1-nht:nc,1-nht:nc)
         !
         !
         ! fill in "e" on Figure above
         !
         ! (use code from west above)
         !
         do halo=1,nhr
           ftmp(:) = fcube(halo,:)   ! copy to a temporary
           do i=halo-nh,min(nc+nh-(halo-1),nc+1)
             ibaseref=ibase(i,halo,1)
             fotherpanel(halo ,i,2) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
           end do
         end do
         !
         !
         ! fill in "n" on Figure above
         !
         ! (use code from north above)
         !
         do halo=1,nhr
           ftmp(:) = fcube(halo,:)   ! copy to a temporary
           do i=0,nht-halo
             ibaseref = ibase(i,halo,2)+nc
             fotherpanel(1-i,nc+halo,2) = dotproduct(hWeight(:,i,halo,2),ftmp(ibaseref:ibaseref+ns-1),ns) !north
           end do
         end do
         fotherpanel(1,nc+1,2)=0.25_r8*(&
              fotherpanel(2,nc+1,2)+fotherpanel(0,nc+1,2)+fotherpanel(1,nc+2,2)+fotherpanel(1,nc,2))
       end if

     else if (cubeboundary==neast) then
       !
       !
       !     0---/---/---/---0---/---/---/-------0     !   x---------------x---------------x--------
       !    0   /   /   /   0   /   /   /      0 |     !   |   |   |   |   | n | n | n | n | n |   |
       !   0---/---/---/---0---/---/---/-----0   |     !   x---------------x---------------x--------
       !  0   /   /   /   0   /   /   /    0 |   |     !   |   |   |   | n | n | n | n | n | a | e |
       ! 000000000000000000000000000000000   |   |     !   x---------------00000000000000000--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 0---------------0---------------0--/|   |     !   x---------------0---------------0--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 0---------------0---------------0--/|   |     !   x---------------0---------------0--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 0---------------0---------------0--/|   0     !   x---------------0---------------0--------
       ! 0   |   |   |   0   |   |   |   0   0000|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 000000000000000000000000000000000000|   |     !   x---------------00000000000000000--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H | r | r | r | r | e | e |
       ! 0---------------0---------------0--/|   |     !   x---------------x---------------x--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   |   | H | H | H | H | H | e |   |
       ! 0---------------0---------------0--/|   |     !   x---------------x---------------x--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   |   |   | H | H | H | H |   |   |
       ! 0---------------0---------------0--/|   0     !   x---------------x---------------x--------
       ! 0   |   |   |   0   |   |   |   0   0000      !   |   |   |   |   |   |   |   |   |   |   |
       ! 000000000000000000000000000000000000          !   x---------------x---------------x--------
       !
       !
       !
       fpanel(1-nht:nc,1-nht:nc)=fcube(1-nht:nc,1-nht:nc)
       !     fotherpanel (nc+1 :nc+nht ,1-nht:nc+nht)=fcube(nc+1 :nc+nht ,1-nht:nc+nht)
       !
       ! east
       !
       do halo=1,nhr
         ftmp(:) = fcube(nc+halo,:)   ! copy to a temporary
         do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=ibase(i,halo,1 )
           fpanel(nc+halo,i) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
         end do
       end do
       !
       ! north
       !
       do halo=1,nhr
         do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=ibase(i,halo,1)
           fpanel(i,nc+halo) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,nc+halo  ),ns) !north
         end do
       end do
       fpanel(nc+1,nc+1)=0.25_r8*(&
            fpanel(nc,nc+1)+fpanel(nc+1,nc)+fpanel(nc+1,nc+2)+fpanel(nc+2,nc+1))
       !
       ! ****************************************************************
       !
       ! fill halo for reconstruction on north neighbor panel projection
       !
       ! ****************************************************************
       !
       ! On the north panel projection the neighbors are arragened as follow (seast case):
       !
       !
       !             |     |
       !             |  N  | E
       !             |-----|
       !                   \
       !                S   \
       !                     \
       !
       !
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   |   |   |   |   |   |   |   |   |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   |   |   | i | i | i | i | i |   |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   |   | i | i | i | i | i | i | e |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   | i | i | r | r | r | r | r | e | e |   |   |
       ! x---x---x---x---00000000000000000---x---x---x---x
       ! |   |   | s | s 0 s | s | s | s 0 se| e |   |   |
       ! x---x---x---x---0---x---x---x---0---x---x---x---x
       ! |   |   |   | s 0 s | s | s | s 0 s |   |   |   |
       ! x---x---x---x---0---x---x---x---0---x---x---x---x
       ! |   |   |   |   0   |   |   |   0   |   |   |   |
       ! x---x---x---x---0---x---x---x---0---x---x---x---x
       ! |   |   |   |   0   |   |   |   0   |   |   |   |
       ! x---x---x---x---00000000000000000---x---x---x---x
       ! |   |   |   |   |   |   |   |   |   |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       !
       !
       ! fill values on same panel projection ("r" and "i" on Figure above)
       !
       if (present(fotherpanel)) then
         fotherpanel(1-nht:nc,nc+1:nc+nht,1)  = fcube(1-nht:nc,nc+1:nc+nht)
         !
         ! fill in "s" on Figure above
         !
         ! (use north case from above and shift/reverse j-index
         !
         do halo=1,nhr
           do i=halo-nh,min(nc+nh-(halo-1),nc+1)
             ibaseref=ibase(i,halo,1)
             fotherpanel (i,nc+1-halo,1) = dotproduct(hWeight(:,i,halo,1),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
           end do
         end do
         !
         ! fill in "e" on Figure above
         !
         do halo=1,nhr
           do i=max(halo-nh,0),nht-halo
             ibaseref=ibase(i,halo,2) +nc
             !
             ! fotherpanel uses indexing of main panel's projection
             ! fcube: rotated indexing
             !
             fotherpanel (nc+halo,nc+i,1) = dotproduct(hWeight(:,i,halo,2),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
           end do
         end do
         fotherpanel(nc+1,nc,1)=0.25_r8*(&
              fotherpanel(nc+2,nc,1)+fotherpanel(nc,nc,1)+fotherpanel(nc+1,nc+1,1)+fotherpanel(nc+1,nc-1,1))
         !
         ! ****************************************************************
         !
         ! fill halo for reconstruction on east neighbor panel projection
         !
         ! ****************************************************************
         !
         ! On the north panel projection the neighbors are arragened as follow (seast case):
         !
         !
         !           \    N
         !            \
         !             \------
         !             |     |
         !         P   |  E  |
         !             |     |
         !             -------
         !
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   | n | n |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   | w | wn| n | n |   |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   |   |   |   |   |   | w | w | r | i | i |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   | w | i | i |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   | i |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---
         !
         !
         !
         ! fill values on same panel projection ("r" and "i" on Figure above)
         !
         fotherpanel(nc+1:nc+nht,1-nht:nc,2)  = fcube(nc+1:nc+nht,1-nht:nc)
         !
         ! fill in "w" on Figure above
         !
         ! (use east case from above and shift/reverse j-index
         !
         do halo=1,nhr
           ftmp(:) = fcube(nc+1-halo,:)   ! copy to a temporary
           do i=halo-nh,min(nc+nh-(halo-1),nc+1)
             ibaseref=ibase(i,halo,1 )
             fotherpanel(nc+1-halo,i,2) = dotproduct(hWeight(:,i,halo,1),ftmp(ibaseref:ibaseref+ns-1),ns)
           end do
         end do
         !
         ! fill in "n" on Figure above
         !
         do halo=1,nhr
           ftmp(:) = fcube(nc+1-halo,:)   ! copy to a temporary
           do i=max(halo-nh,0),nht-halo
             ibaseref=ibase(i,halo,2) +nc
             !
             ! fotherpanel uses indexing of main panel's projection
             ! fcube: rotated indexing
             !
             fotherpanel (nc+i,nc+halo,2) = dotproduct(hWeight(:,i,halo,2),ftmp(ibaseref:ibaseref+ns-1),ns)
           end do
         end do
         fotherpanel(nc,nc+1,2)=0.25_r8*(&
              fotherpanel(nc+1,nc+1,2)+fotherpanel(nc-1,nc+1,2)+fotherpanel(nc,nc+2,2)+fotherpanel(nc,nc,2))
       end if
     end if
   end subroutine extend_panel_interpolate
   !
   ! initialize non-existent ghost cells
   !
   subroutine zero_non_existent_ghost_cell(recons,irecons,cubeboundary,nc,nhe,ntrac_in)
     use control_mod, only : north, south, east, west, neast, nwest, seast, swest

     integer,          intent(in)  :: nc,nhe,cubeboundary,irecons,ntrac_in
     real (kind=r8), dimension(irecons,1-nhe:nc+nhe,1-nhe:nc+nhe,ntrac_in), intent(out):: recons

     integer :: i,j

     if (cubeboundary>0) then
        if (cubeboundary==nwest) then
           do j=nc+1,nc+nhe
              do i=1-nhe,0
                 recons(:,i,j,:) = 0.0_r8
              end do
           end do
        else if (cubeboundary==swest) then
           do j=1-nhe,0
              do i=1-nhe,0
                 recons(:,i,j,:) = 0.0_r8
              end do
           end do
        else if (cubeboundary==neast) then
           do j=nc+1,nc+nhe
              do i=nc+1,nc+nhe
                 recons(:,i,j,:) = 0.0_r8
              end do
           end do
        else if (cubeboundary==seast) then
           do j=1-nhe,0
              do i=nc+1,nc+nhe
                 recons(:,i,j,:) = 0.0_r8
              end do
           end do
        end if
     end if
   end subroutine zero_non_existent_ghost_cell
end module fvm_reconstruction_mod
