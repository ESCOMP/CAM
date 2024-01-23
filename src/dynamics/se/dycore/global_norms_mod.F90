!#define old_cam
module global_norms_mod

  use shr_kind_mod, only: r8=>shr_kind_r8
  use cam_logfile,  only: iulog
  use edgetype_mod, only: EdgeBuffer_t

  implicit none
  private
  save

  public :: l1_snorm
  public :: l2_snorm
  public :: linf_snorm

  public :: l1_vnorm
  public :: l2_vnorm
  public :: linf_vnorm

  public :: print_cfl
  public :: global_integral
  public :: global_integrals_general
  public :: wrap_repro_sum

  private :: global_maximum
  type (EdgeBuffer_t), private :: edgebuf

  interface global_integral
     module procedure global_integral_elem
     module procedure global_integral_fvm
  end interface global_integral

contains


  subroutine global_integrals(elem,fld,hybrid,npts,num_flds,nets,nete,I_sphere)
    use hybrid_mod,     only: hybrid_t
    use element_mod,    only: element_t
    use dimensions_mod, only: np
    use physconst,      only: pi
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete,num_flds
    real (kind=r8), intent(in) :: fld(npts,npts,num_flds,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=r8) :: I_sphere(num_flds)
    !
    ! Local variables
    !
    integer :: ie,j,i,q

    real (kind=r8) :: da
    real (kind=r8) :: J_tmp(nets:nete,num_flds)
    !
    ! This algorithm is independent of thread count and task count.
    ! This is a requirement of consistancy checking in cam.
    !
    J_tmp = 0.0_r8

    do ie=nets,nete
      do q=1,num_flds
        do j=1,np
          do i=1,np
            da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
            J_tmp(ie,q) = J_tmp(ie,q) + da*fld(i,j,q,ie)
          end do
        end do
      end do
    end do
    do ie=nets,nete
      global_shared_buf(ie,1:num_flds) = J_tmp(ie,:)
    enddo
    call wrap_repro_sum(nvars=num_flds, comm=hybrid%par%comm)
    I_sphere(:) =global_shared_sum(1:num_flds) /(4.0_r8*PI)
  end subroutine global_integrals

  subroutine global_integrals_general(fld,hybrid,npts,da,num_flds,nets,nete,I_sphere)
    use hybrid_mod,     only: hybrid_t
    use physconst,      only: pi
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    integer,         intent(in) :: npts,nets,nete,num_flds
    real (kind=r8),  intent(in) :: fld(npts,npts,num_flds,nets:nete)
    type (hybrid_t), intent(in) :: hybrid
    real (kind=r8),  intent(in) :: da(npts,npts,nets:nete)

    real (kind=r8) :: I_sphere(num_flds)
    !
    ! Local variables
    !
    integer :: ie,j,i,q

    real (kind=r8) :: J_tmp(nets:nete,num_flds)
    !
    ! This algorithm is independent of thread count and task count.
    ! This is a requirement of consistancy checking in cam.
    !
    J_tmp = 0.0_r8

    do ie=nets,nete
      do q=1,num_flds
        do j=1,npts
          do i=1,npts
            J_tmp(ie,q) = J_tmp(ie,q) + da(i,j,ie)*fld(i,j,q,ie)
          end do
        end do
      end do
    end do
    do ie=nets,nete
      global_shared_buf(ie,1:num_flds) = J_tmp(ie,:)
    enddo
    call wrap_repro_sum(nvars=num_flds, comm=hybrid%par%comm)
    I_sphere(:) =global_shared_sum(1:num_flds) /(4.0_r8*PI)
  end subroutine global_integrals_general


  ! ================================
  ! global_integral:
  !
  ! eq 81 in Williamson, et. al. p 218
  ! for spectral elements
  !
  ! ================================
  ! --------------------------
  function global_integral_elem(elem,fld,hybrid,npts,nets,nete) result(I_sphere)
    use hybrid_mod,     only: hybrid_t
    use element_mod,    only: element_t
    use dimensions_mod, only: np
    use physconst,      only: pi
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: fld(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=r8) :: I_sphere

    ! Local variables

    integer :: ie,j,i
    real(kind=r8) :: I_tmp(1)

    real (kind=r8) :: da
    real (kind=r8) :: J_tmp(nets:nete)
!
! This algorithm is independent of thread count and task count.
! This is a requirement of consistancy checking in cam.
!
    J_tmp = 0.0_r8

       do ie=nets,nete
          do j=1,np
             do i=1,np
                da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
                J_tmp(ie) = J_tmp(ie) + da*fld(i,j,ie)
             end do
          end do
       end do
    do ie=nets,nete
      global_shared_buf(ie,1) = J_tmp(ie)
    enddo
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
    I_tmp = global_shared_sum(1)
    I_sphere = I_tmp(1)/(4.0_r8*PI)

  end function global_integral_elem

  function global_integral_fvm(fvm,fld,hybrid,npts,nets,nete) result(I_sphere)
    use hybrid_mod,     only: hybrid_t
    use fvm_control_volume_mod, only: fvm_struct
    use physconst,      only: pi
    use parallel_mod,   only: global_shared_buf, global_shared_sum

    type (fvm_struct)    , intent(in) :: fvm(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: fld(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=r8) :: I_sphere

    ! Local variables

    integer :: ie,j,i
    real(kind=r8) :: I_tmp(1)

    real (kind=r8) :: da
    real (kind=r8) :: J_tmp(nets:nete)
!
! This algorithm is independent of thread count and task count.
! This is a requirement of consistancy checking in cam.
!
    J_tmp = 0.0_r8
    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             da = fvm(ie)%area_sphere(i,j)
             J_tmp(ie) = J_tmp(ie) + da*fld(i,j,ie)
          end do
       end do
    end do
    do ie=nets,nete
       global_shared_buf(ie,1) = J_tmp(ie)
    enddo
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
    I_tmp = global_shared_sum(1)
    I_sphere = I_tmp(1)/(4.0_r8*PI)

  end function global_integral_fvm

!------------------------------------------------------------------------------------

  ! ================================
  ! print_cfl:
  !
  ! Calculate / output CFL info
  ! (both advective and based on
  ! viscosity or hyperviscosity)
  !
  ! ================================

  subroutine print_cfl(elem,hybrid,nets,nete,dtnu,ptop,pmid,&
       dt_remap_actual,dt_tracer_fvm_actual,dt_tracer_se_actual,&
       dt_dyn_actual,dt_dyn_visco_actual,dt_dyn_del2_actual,dt_tracer_visco_actual,dt_phys)
    !
    !   estimate various CFL limits
    !   also, for variable resolution viscosity coefficient, make sure
    !   worse viscosity CFL (given by dtnu) is not violated by reducing
    !   viscosity coefficient in regions where CFL is violated
    !
    use hybrid_mod,     only: hybrid_t
    use element_mod,    only: element_t
    use dimensions_mod, only: np,ne,nelem,nc,nhe,use_cslam,nlev,large_Courant_incr
    use dimensions_mod, only: nu_scale_top,nu_div_lev,nu_lev,nu_t_lev

    use quadrature_mod, only: gausslobatto, quadrature_t

    use reduction_mod,  only: ParallelMin,ParallelMax
    use physconst,      only: ra, rearth, pi
    use control_mod,    only: nu, nu_div, nu_q, nu_p, nu_t, nu_top, fine_ne, max_hypervis_courant
    use control_mod,    only: tstep_type, hypervis_power, hypervis_scaling
    use control_mod,    only: sponge_del4_nu_div_fac, sponge_del4_nu_fac, sponge_del4_lev
    use cam_abortutils, only: endrun
    use parallel_mod,   only: global_shared_buf, global_shared_sum
    use edge_mod,       only: initedgebuffer, FreeEdgeBuffer, edgeVpack, edgeVunpack
    use bndry_mod,      only: bndry_exchange
    use mesh_mod,       only: MeshUseMeshFile
    use dimensions_mod, only: ksponge_end, kmvis_ref, kmcnd_ref,rho_ref
    use physconst,      only: cpair
    use std_atm_profile,only: std_atm_height
    type(element_t)      , intent(inout) :: elem(:)
    integer              , intent(in) :: nets,nete
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8), intent(in) :: dtnu, ptop, pmid(nlev)
    !
    ! actual time-steps
    !
    real (kind=r8), intent(in) :: dt_remap_actual,dt_tracer_fvm_actual,dt_tracer_se_actual,&
                           dt_dyn_actual,dt_dyn_visco_actual,dt_dyn_del2_actual,           &
                           dt_tracer_visco_actual, dt_phys

    ! Element statisics
    real (kind=r8) :: max_min_dx,min_min_dx,min_max_dx,max_unif_dx   ! used for normalizing scalar HV
    real (kind=r8) :: max_normDinv, min_normDinv  ! used for CFL
    real (kind=r8) :: min_area, max_area,max_ratio !min/max element area
    real (kind=r8) :: avg_area, avg_min_dx,tot_area,tot_area_rad
    real (kind=r8) :: min_hypervis, max_hypervis, avg_hypervis, stable_hv
    real (kind=r8) :: normDinv_hypervis
    real (kind=r8) :: x, y, noreast, nw, se, sw
    real (kind=r8), dimension(np,np,nets:nete) :: zeta
    real (kind=r8) :: lambda_max, lambda_vis, min_gw, lambda,umax, ugw
    real (kind=r8) :: scale1, max_laplace,z(nlev)
    integer :: ie, i, j, rowind, colind, k
    type (quadrature_t)    :: gp
    character(LEN=256) :: rk_str

    real (kind=r8) :: s_laplacian, s_hypervis, s_rk, s_rk_tracer !Stability region
    real (kind=r8) :: dt_max_adv, dt_max_gw, dt_max_tracer_se, dt_max_tracer_fvm
    real (kind=r8) :: dt_max_hypervis, dt_max_hypervis_tracer, dt_max_laplacian_top

    real(kind=r8) :: I_sphere, nu_max, nu_div_max
    real(kind=r8) :: fld(np,np,nets:nete)

    logical :: top_000_032km, top_032_042km, top_042_090km, top_090_140km, top_140_600km ! model top location ranges
    logical :: nu_set,div_set,lev_set

    ! Eigenvalues calculated by folks at UMich (Paul U & Jared W)
    select case (np)
    case (2)
      lambda_max = 0.5_r8
      lambda_vis = 0.0_r8  ! need to compute this
    case (3)
      lambda_max = 1.5_r8
      lambda_vis = 12.0_r8
    case (4)
      lambda_max = 2.74_r8
      lambda_vis = 30.0_r8
    case (5)
      lambda_max = 4.18_r8
      lambda_vis = 91.6742_r8
    case (6)
      lambda_max = 5.86_r8
      lambda_vis = 190.1176_r8
    case (7)
      lambda_max = 7.79_r8
      lambda_vis = 374.7788_r8
    case (8)
      lambda_max = 10.0_r8
      lambda_vis = 652.3015_r8
    case DEFAULT
      lambda_max = 0.0_r8
      lambda_vis = 0.0_r8
    end select

    if ((lambda_max.eq.0_r8).and.(hybrid%masterthread)) then
      print*, "lambda_max not calculated for NP = ",np
      print*, "Estimate of gravity wave timestep will be incorrect"
    end if
    if ((lambda_vis.eq.0_r8).and.(hybrid%masterthread)) then
      print*, "lambda_vis not calculated for NP = ",np
      print*, "Estimate of viscous CFLs will be incorrect"
    end if

    do ie=nets,nete
      elem(ie)%variable_hyperviscosity = 1.0_r8
    end do

    gp=gausslobatto(np)
    min_gw = minval(gp%weights)
    !
    !******************************************************************************************
    !
    ! compute some local and global grid metrics
    !
    !******************************************************************************************
    !
    fld(:,:,nets:nete)=1.0_r8
    ! Calculate surface area by integrating 1.0_r8 over sphere and dividing by 4*PI (Should be 1)
    I_sphere = global_integral(elem, fld(:,:,nets:nete),hybrid,np,nets,nete)

    min_normDinv = 1E99_r8
    max_normDinv = 0
    min_max_dx   = 1E99_r8
    min_min_dx   = 1E99_r8
    max_min_dx   = 0
    min_area     = 1E99_r8
    max_area     = 0
    max_ratio    = 0
    do ie=nets,nete
      max_normDinv  = max(max_normDinv,elem(ie)%normDinv)
      min_normDinv  = min(min_normDinv,elem(ie)%normDinv)
      min_min_dx    = min(min_min_dx,elem(ie)%dx_short)
      max_min_dx    = max(max_min_dx,elem(ie)%dx_short)
      min_max_dx    = min(min_max_dx,elem(ie)%dx_long)

      elem(ie)%area = sum(elem(ie)%spheremp(:,:))
      min_area      = min(min_area,elem(ie)%area)
      max_area      = max(max_area,elem(ie)%area)
      max_ratio     = max(max_ratio,elem(ie)%dx_long/elem(ie)%dx_short)

      global_shared_buf(ie,1) = elem(ie)%area
      global_shared_buf(ie,2) = elem(ie)%dx_short
    enddo
    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)
    avg_area     = global_shared_sum(1)/dble(nelem)
    tot_area_rad = global_shared_sum(1)
    avg_min_dx   = global_shared_sum(2)/dble(nelem)

    min_area     = ParallelMin(min_area,hybrid)
    max_area     = ParallelMax(max_area,hybrid)
    min_normDinv = ParallelMin(min_normDinv,hybrid)
    max_normDinv = ParallelMax(max_normDinv,hybrid)
    min_min_dx   = ParallelMin(min_min_dx,hybrid)
    max_min_dx   = ParallelMax(max_min_dx,hybrid)
    min_max_dx   = ParallelMin(min_max_dx,hybrid)
    max_ratio    = ParallelMax(max_ratio,hybrid)
    ! Physical units for area (unit sphere to Earth sphere)
    min_area = min_area*rearth*rearth/1000000._r8    !m2 (rearth is in units of km)
    max_area = max_area*rearth*rearth/1000000._r8    !m2 (rearth is in units of km)
    avg_area = avg_area*rearth*rearth/1000000._r8    !m2 (rearth is in units of km)
    tot_area = tot_area_rad*rearth*rearth/1000000._r8!m2 (rearth is in units of km)
    if (hybrid%masterthread) then
       write(iulog,* )""
       write(iulog,* )"Running Global Integral Diagnostic..."
       write(iulog,*)"Area of unit sphere is",I_sphere
       write(iulog,*)"Should be 1.0 to round off..."
       write(iulog,'(a,f9.3)') 'Element area:  max/min',(max_area/min_area)
       write(iulog,'(a,E23.15)') 'Total Grid area:  ',(tot_area)
       write(iulog,'(a,E23.15)') 'Total Grid area rad^2:  ',(tot_area_rad)
       if (.not.MeshUseMeshFile) then
           write(iulog,'(a,f6.3,f8.2)') "Average equatorial node spacing (deg, km) = ", &
                dble(90)/dble(ne*(np-1)), PI*rearth/(2000.0_r8*dble(ne*(np-1)))
       end if
       write(iulog,'(a,2f9.3)') 'norm of Dinv (min, max): ', min_normDinv, max_normDinv
       write(iulog,'(a,1f8.2)') 'Max Dinv-based element distortion: ', max_ratio
       write(iulog,'(a,3f8.2)') 'dx based on Dinv svd:          ave,min,max = ', avg_min_dx, min_min_dx, max_min_dx
       write(iulog,'(a,3f8.2)') "dx based on sqrt element area: ave,min,max = ", &
                sqrt(avg_area)/(np-1),sqrt(min_area)/(np-1),sqrt(max_area)/(np-1)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  SCALAR, RESOLUTION-AWARE HYPERVISCOSITY
    !  this block of code initializes the variable_hyperviscsoity() array
    !  based on largest length scale in each element and user specified scaling
    !  it then limits the coefficient if the user specifed a max CFL
    !  this limiting is based on the smallest length scale of each element
    !  since that controls the CFL.
    !  Mike Levy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (hypervis_power /= 0) then

      min_hypervis = 1d99
      max_hypervis = 0
      avg_hypervis = 0


      max_unif_dx = min_max_dx  ! use this for average resolution, unless:
      !       viscosity in namelist specified for smallest element:
      if (fine_ne>0) then
        ! viscosity in namelist specified for regions with a resolution
        ! equivilant to a uniform grid with ne=fine_ne
        if (np /= 4 ) call endrun('ERROR: setting fine_ne only supported with NP=4')
        max_unif_dx = (111.28_r8*30)/dble(fine_ne)   ! in km
      endif

      !
      ! note: if L = eigenvalue of metinv, then associated length scale (km) is
      ! dx = 1.0_r8/( sqrt(L)*0.5_r8*dble(np-1)*ra*1000.0_r8)
      !
      !       for viscosity *tensor*, we take at each point:
      !            nu1 = nu*(dx1/max_unif_dx)**3.2      dx1 associated with eigenvalue 1
      !            nu2 = nu*(dx2/max_unif_dx)**3.2      dx2 associated with eigenvalue 2
      !       with this approach:
      !          - with this formula, no need to adjust for CFL violations
      !          - if nu comes from a 3.2 scaling that is stable for coarse and fine resolutions,
      !            this formulat will be stable.
      !          - gives the correct answer in long skinny rectangles:
      !            large viscosity in the long direction, small viscosity in the short direction
      !
      normDinv_hypervis = 0
      do ie=nets,nete
        ! variable viscosity based on map from ulatlon -> ucontra

        ! dx_long
        elem(ie)%variable_hyperviscosity = sqrt((elem(ie)%dx_long/max_unif_dx) ** hypervis_power)
        elem(ie)%hv_courant = dtnu*(elem(ie)%variable_hyperviscosity(1,1)**2) * &
             (lambda_vis**2) * ((ra*elem(ie)%normDinv)**4)

        ! Check to see if this is stable
        if (elem(ie)%hv_courant.gt.max_hypervis_courant) then
          stable_hv = sqrt( max_hypervis_courant / &
               (  dtnu * (lambda_vis)**2 * (ra*elem(ie)%normDinv)**4 ) )

#if 0
          ! Useful print statements for debugging the adjustments to hypervis
          print*, "Adjusting hypervis on elem ", elem(ie)%GlobalId
          print*, "From ", nu*elem(ie)%variable_hyperviscosity(1,1)**2, " to ", nu*stable_hv
          print*, "Difference = ", nu*(/elem(ie)%variable_hyperviscosity(1,1)**2-stable_hv/)
          print*, "Factor of ", elem(ie)%variable_hyperviscosity(1,1)**2/stable_hv
          print*, " "
#endif
          !                make sure that: elem(ie)%hv_courant <=  max_hypervis_courant
          elem(ie)%variable_hyperviscosity = stable_hv
          elem(ie)%hv_courant = dtnu*(stable_hv**2) * (lambda_vis)**2 * (ra*elem(ie)%normDinv)**4
        end if
        normDinv_hypervis = max(normDinv_hypervis, elem(ie)%hv_courant/dtnu)

        min_hypervis = min(min_hypervis, elem(ie)%variable_hyperviscosity(1,1))
        max_hypervis = max(max_hypervis, elem(ie)%variable_hyperviscosity(1,1))
        global_shared_buf(ie,1) = elem(ie)%variable_hyperviscosity(1,1)
      end do

      min_hypervis = ParallelMin(min_hypervis, hybrid)
      max_hypervis = ParallelMax(max_hypervis, hybrid)
      call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
      avg_hypervis = global_shared_sum(1)/dble(nelem)

      normDinv_hypervis = ParallelMax(normDinv_hypervis, hybrid)

      ! apply DSS (aka assembly procedure) to variable_hyperviscosity (makes continuous)
      call initEdgeBuffer(hybrid%par,edgebuf,elem,1)
      do ie=nets,nete
        zeta(:,:,ie) = elem(ie)%variable_hyperviscosity(:,:)*elem(ie)%spheremp(:,:)
        call edgeVpack(edgebuf,zeta(1,1,ie),1,0,ie)
      end do
      call bndry_exchange(hybrid,edgebuf,location='print_cfl #1')
      do ie=nets,nete
        call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,ie)
        elem(ie)%variable_hyperviscosity(:,:) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
      end do
      call FreeEdgeBuffer(edgebuf)

      ! replace hypervis w/ bilinear based on continuous corner values
      do ie=nets,nete
        noreast = elem(ie)%variable_hyperviscosity(np,np)
        nw = elem(ie)%variable_hyperviscosity(1,np)
        se = elem(ie)%variable_hyperviscosity(np,1)
        sw = elem(ie)%variable_hyperviscosity(1,1)
        do i=1,np
          x = gp%points(i)
          do j=1,np
            y = gp%points(j)
            elem(ie)%variable_hyperviscosity(i,j) = 0.25_r8*( &
                 (1.0_r8-x)*(1.0_r8-y)*sw + &
                 (1.0_r8-x)*(y+1.0_r8)*nw + &
                 (x+1.0_r8)*(1.0_r8-y)*se + &
                 (x+1.0_r8)*(y+1.0_r8)*noreast)
          end do
        end do
      end do
    else  if (hypervis_scaling/=0) then
      ! tensorHV.  New eigenvalues are the eigenvalues of the tensor V
      ! formulas here must match what is in cube_mod.F90
      ! for tensorHV, we scale out the rearth dependency
      lambda = max_normDinv**2
      normDinv_hypervis = (lambda_vis**2) * (max_normDinv**4) * &
           (lambda**(-hypervis_scaling/2) )
    else
      ! constant coefficient formula:
      normDinv_hypervis = (lambda_vis**2) * (ra*max_normDinv)**4
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  TENSOR, RESOLUTION-AWARE HYPERVISCOSITY
    !  The tensorVisc() array is computed in cube_mod.F90
    !  this block of code will DSS it so the tensor if C0
    !  and also make it bilinear in each element.
    !  Oksana Guba
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (hypervis_scaling /= 0) then

      call initEdgeBuffer(hybrid%par,edgebuf,elem,1)
      do rowind=1,2
        do colind=1,2
          do ie=nets,nete
            zeta(:,:,ie) = elem(ie)%tensorVisc(:,:,rowind,colind)*elem(ie)%spheremp(:,:)
            call edgeVpack(edgebuf,zeta(1,1,ie),1,0,ie)
          end do

          call bndry_exchange(hybrid,edgebuf)
          do ie=nets,nete
            call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,ie)
            elem(ie)%tensorVisc(:,:,rowind,colind) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
          end do
        enddo !rowind
      enddo !colind
      call FreeEdgeBuffer(edgebuf)

      !IF BILINEAR MAP OF V NEEDED

      do rowind=1,2
        do colind=1,2
          ! replace hypervis w/ bilinear based on continuous corner values
          do ie=nets,nete
            noreast = elem(ie)%tensorVisc(np,np,rowind,colind)
            nw = elem(ie)%tensorVisc(1,np,rowind,colind)
            se = elem(ie)%tensorVisc(np,1,rowind,colind)
            sw = elem(ie)%tensorVisc(1,1,rowind,colind)
            do i=1,np
              x = gp%points(i)
              do j=1,np
                y = gp%points(j)
                elem(ie)%tensorVisc(i,j,rowind,colind) = 0.25_r8*( &
                     (1.0_r8-x)*(1.0_r8-y)*sw + &
                     (1.0_r8-x)*(y+1.0_r8)*nw + &
                     (x+1.0_r8)*(1.0_r8-y)*se + &
                     (x+1.0_r8)*(y+1.0_r8)*noreast)
              end do
            end do
          end do
        enddo !rowind
      enddo !colind
    endif
    deallocate(gp%points)
    deallocate(gp%weights)

    call automatically_set_viscosity_coefficients(hybrid,ne,max_min_dx,min_min_dx,nu_p  ,1.0_r8 ,'_p  ')
    call automatically_set_viscosity_coefficients(hybrid,ne,max_min_dx,min_min_dx,nu    ,1.0_r8,'    ')
    call automatically_set_viscosity_coefficients(hybrid,ne,max_min_dx,min_min_dx,nu_div,2.5_r8 ,'_div')

    if (nu_q<0) nu_q = nu_p ! necessary for consistency
    if (nu_t<0) nu_t = nu_p ! temperature damping is always equal to nu_p

    nu_div_lev(:) = nu_div
    nu_lev(:)     = nu
    nu_t_lev(:)   = nu_p

    !
    ! sponge layer strength needed for stability depends on model top location
    !
    top_000_032km = .false.
    top_032_042km = .false.
    top_042_090km = .false.
    top_090_140km = .false.
    top_140_600km = .false.
    nu_set = sponge_del4_nu_fac < 0
    div_set = sponge_del4_nu_div_fac < 0
    lev_set = sponge_del4_lev < 0
    if (ptop>1000.0_r8) then
      !
      ! low top (~10 Pa)
      !
      top_000_032km = .true.
    else if (ptop>100.0_r8) then
      !
      ! CAM6 top (~225 Pa) or CAM7 low top
      !
      top_032_042km = .true.
    else if (ptop>1e-1_r8) then
      !
      ! CAM7 top (~4.35e-1 Pa)
      !
      top_042_090km = .true.
    else if (ptop>1E-4_r8) then
      !
      ! WACCM top (~4.5e-4 Pa)
      !
      top_090_140km = .true.
    else
      !
      ! WACCM-x - geospace (~4e-7 Pa)
      !
      top_140_600km = .true.
    end if
    !
    ! Logging text for sponge layer configuration
    !
    if (hybrid%masterthread .and. (nu_set .or. div_set .or. lev_set)) then
       write(iulog,* )""
       write(iulog,* )"Sponge layer del4 coefficient defaults based on model top location:"
    end if
    !
    ! if user or namelist is not specifying sponge del4 settings here are best guesses (empirically determined)
    !
    umax = 0.0_r8
    if (top_000_032km) then
      umax = 120._r8
      if (sponge_del4_lev       <0) sponge_del4_lev        = 1
      if (sponge_del4_nu_fac    <0) sponge_del4_nu_fac     = 1.0_r8
      if (sponge_del4_nu_div_fac<0) sponge_del4_nu_div_fac = 1.0_r8
    end if

    if (top_032_042km) then
      umax = 120._r8
      if (sponge_del4_lev       <0) sponge_del4_lev        = 1
      if (sponge_del4_nu_fac    <0) sponge_del4_nu_fac     = 1.0_r8
      if (sponge_del4_nu_div_fac<0) sponge_del4_nu_div_fac = 1.0_r8
    end if

    if (top_042_090km) then
      umax = 240._r8
      if (sponge_del4_lev       <0) sponge_del4_lev        = 1
      if (sponge_del4_nu_fac    <0) sponge_del4_nu_fac     = 1.0_r8
      if (sponge_del4_nu_div_fac<0) sponge_del4_nu_div_fac = 1.0_r8
    end if

    if (top_090_140km) then
      umax = 300._r8
    end if
    if (top_140_600km) then
      umax = 800._r8
    end if
    if (top_090_140km.or.top_140_600km) then
      if (sponge_del4_lev       <0) sponge_del4_lev        = 20
      if (sponge_del4_nu_fac    <0) sponge_del4_nu_fac     = 5.0_r8
      if (sponge_del4_nu_div_fac<0) sponge_del4_nu_div_fac = 10.0_r8
    end if
    !
    ! Log sponge layer configuration
    !
    if (hybrid%masterthread) then
       if (nu_set) then
          write(iulog, '(a,e9.2)')   '  sponge_del4_nu_fac     = ',sponge_del4_nu_fac
       end if

       if (div_set) then
          write(iulog, '(a,e9.2)')   '  sponge_del4_nu_div_fac = ',sponge_del4_nu_div_fac
       end if

       if (lev_set) then
          write(iulog, '(a,i0)')   '  sponge_del4_lev        = ',sponge_del4_lev
       end if
       write(iulog,* )""
    end if

    nu_max     =  sponge_del4_nu_fac*nu_p
    nu_div_max =  sponge_del4_nu_div_fac*nu_p
    do k=1,nlev
      ! Vertical profile from FV dycore (see Lauritzen et al. 2012 DOI:10.1177/1094342011410088)
      scale1        = 0.5_r8*(1.0_r8+tanh(2.0_r8*log(pmid(sponge_del4_lev)/pmid(k))))
      if (sponge_del4_nu_div_fac /= 1.0_r8) then
        nu_div_lev(k) = (1.0_r8-scale1)*nu_div+scale1*nu_div_max
      end if
      if (sponge_del4_nu_fac /= 1.0_r8) then
        nu_lev(k)     = (1.0_r8-scale1)*nu    +scale1*nu_max
        nu_t_lev(k)   = (1.0_r8-scale1)*nu_p  +scale1*nu_max
      end if
    end do

    if (hybrid%masterthread)then
      write(iulog,*) "z computed from barometric formula (using US std atmosphere)"
      call std_atm_height(pmid(:),z(:))
      write(iulog,*) "k,pmid_ref,z,nu_lev,nu_t_lev,nu_div_lev"
      do k=1,nlev
        write(iulog,'(i3,5e11.4)') k,pmid(k),z(k),nu_lev(k),nu_t_lev(k),nu_div_lev(k)
      end do
      if (nu_top>0) then
        write(iulog,*) ": ksponge_end = ",ksponge_end
        write(iulog,*) ": sponge layer Laplacian damping"
        write(iulog,*) "k, p, z, nu_scale_top, nu (actual Laplacian damping coefficient)"

        do k=1,ksponge_end
           write(iulog,'(i3,4e11.4)') k,pmid(k),z(k),nu_scale_top(k),nu_scale_top(k)*nu_top
        end do
      end if
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! time-step information
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! S=time-step stability region (i.e. advection w/leapfrog: S=1, viscosity w/forward Euler: S=2)
    !
    if (tstep_type==1) then
      S_rk   = 2.0_r8
      rk_str = '  * RK2-SSP 3 stage (same as tracers)'
    elseif (tstep_type==2) then
      S_rk   = 2.0_r8
      rk_str = '  * classic RK3'
    elseif (tstep_type==3) then
      S_rk   = 2.0_r8
      rk_str = '  * Kinnmark&Gray RK4'
    elseif (tstep_type==4) then
      S_rk   = 3.0_r8
      rk_str = '  * Kinnmark&Gray RK3 5 stage (3rd order)'
    end if
    if (hybrid%masterthread) then
      write(iulog,'(a,f12.8,a)') 'Model top is ',ptop,'Pa'
      write(iulog,'(a)') ' '
      write(iulog,'(a)') 'Timestepping methods used in dynamical core:'
      write(iulog,'(a)')
      write(iulog,*) rk_str
      write(iulog,'(a)') '   * Spectral-element advection uses SSP preservation RK3'
      write(iulog,'(a)') '   * Viscosity operators use forward Euler'
    end if
    S_laplacian = 2.0_r8 !using forward Euler for sponge diffusion
    S_hypervis  = 2.0_r8 !using forward Euler for hyperviscosity
    S_rk_tracer = 2.0_r8

    ugw = 342.0_r8 !max gravity wave speed

    dt_max_adv             = S_rk/(umax*max_normDinv*lambda_max*ra)
    dt_max_gw              = S_rk/(ugw*max_normDinv*lambda_max*ra)
    dt_max_tracer_se       = S_rk_tracer*min_gw/(umax*max_normDinv*ra)
    if (use_cslam) then
      if (large_Courant_incr) then
        dt_max_tracer_fvm      = dble(nhe)*(4.0_r8*pi*Rearth/dble(4.0_r8*ne*nc))/umax
      else
        dt_max_tracer_fvm      = dble(nhe)*(2.0_r8*pi*Rearth/dble(4.0_r8*ne*nc))/umax
      end if
    else
      dt_max_tracer_fvm = -1.0_r8
    end if
    nu_max = MAX(MAXVAL(nu_div_lev(:)),MAXVAL(nu_lev(:)),MAXVAL(nu_t_lev(:)))
    dt_max_hypervis        = s_hypervis/(nu_max*normDinv_hypervis)
    dt_max_hypervis_tracer = s_hypervis/(nu_q*normDinv_hypervis)

    max_laplace = MAX(MAXVAL(nu_scale_top(:))*nu_top,MAXVAL(kmvis_ref(:)/rho_ref(:)))
    max_laplace = MAX(max_laplace,MAXVAL(kmcnd_ref(:)/(cpair*rho_ref(:))))
    dt_max_laplacian_top   = 1.0_r8/(max_laplace*((ra*max_normDinv)**2)*lambda_vis)

    if (hybrid%masterthread) then
      write(iulog,'(a,f10.2,a)') ' '
      write(iulog,'(a,f10.2,a)') 'Estimates for maximum stable and actual time-steps for different aspects of algorithm:'
      write(iulog,'(a,f12.8,a)') '(assume max wind is ',umax,'m/s)'
      write(iulog,'(a)')         '(assume max gravity wave speed is 342m/s)'
      write(iulog,'(a,f10.2,a)') ' '
      write(iulog,'(a,f10.2,a,f10.2,a)') '* dt_dyn        (time-stepping dycore  ; u,v,T,dM) < ',&
           MIN(dt_max_adv,dt_max_gw),'s ',dt_dyn_actual,'s'
      if (dt_dyn_actual>MIN(dt_max_adv,dt_max_gw)) write(iulog,*) 'WARNING: dt_dyn theoretically unstable'

      write(iulog,'(a,f10.2,a,f10.2,a)') '* dt_dyn_vis    (hyperviscosity)       ; u,v,T,dM) < ',dt_max_hypervis,&
           's ',dt_dyn_visco_actual,'s'
      if (dt_dyn_visco_actual>dt_max_hypervis) write(iulog,*) 'WARNING: dt_dyn_vis theoretically unstable'
      if (.not.use_cslam) then
         write(iulog,'(a,f10.2,a,f10.2,a)') '* dt_tracer_se  (time-stepping tracers ; q       ) < ',dt_max_tracer_se,'s ',&
           dt_tracer_se_actual,'s'
         if (dt_tracer_se_actual>dt_max_tracer_se) write(iulog,*) 'WARNING: dt_tracer_se theoretically unstable'
         write(iulog,'(a,f10.2,a,f10.2,a)') '* dt_tracer_vis (hyperviscosity tracers; q       ) < ',dt_max_hypervis_tracer,'s',&
              dt_tracer_visco_actual,'s'
         if (dt_tracer_visco_actual>dt_max_hypervis_tracer) write(iulog,*) 'WARNING: dt_tracer_hypervis theoretically unstable'
      end if
      if (use_cslam) then
        write(iulog,'(a,f10.2,a,f10.2,a)') '* dt_tracer_fvm (time-stepping tracers ; q       ) < ',dt_max_tracer_fvm,&
             's ',dt_tracer_fvm_actual
        if (dt_tracer_fvm_actual>dt_max_tracer_fvm) write(iulog,*) 'WARNING: dt_tracer_fvm theortically unstable'
      end if
      write(iulog,'(a,f10.2)') '* dt_remap (vertical remap dt) ',dt_remap_actual
      do k=1,ksponge_end
        max_laplace = MAX(nu_scale_top(k)*nu_top,kmvis_ref(k)/rho_ref(k))
        max_laplace = MAX(max_laplace,kmcnd_ref(k)/(cpair*rho_ref(k)))
        dt_max_laplacian_top   = 1.0_r8/(max_laplace*((ra*max_normDinv)**2)*lambda_vis)

        write(iulog,'(a,f10.2,a,f10.2,a)') '* dt    (del2 sponge           ; u,v,T,dM) < ',&
             dt_max_laplacian_top,'s',dt_dyn_del2_actual,'s'
        if (dt_dyn_del2_actual>dt_max_laplacian_top) then
          if (k==1) then
            write(iulog,*) 'WARNING: theoretically unstable in sponge; increase se_hypervis_subcycle_sponge',&
                           ' (this WARNING can sometimes be ignored in level 1)'
          else
            write(iulog,*) 'WARNING: theoretically unstable in sponge; increase se_hypervis_subcycle_sponge'
          endif
        end if
      end do
      write(iulog,*) ' '
      if (hypervis_power /= 0) then
        write(iulog,'(a,3e11.4)')'Scalar hyperviscosity (dynamics): ave,min,max = ', &
             nu*(/avg_hypervis**2,min_hypervis**2,max_hypervis**2/)
      end if
      write(iulog,*) 'tstep_type = ',tstep_type
    end if
  end subroutine print_cfl

  !
  ! ============================
  ! global_maximum:
  !
  ! Find global maximum on sphere
  !
  ! ================================

  function global_maximum(fld,hybrid,npts,nets,nete) result(Max_sphere)

    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : red_max, pmax_mt

    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: fld(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=r8) :: Max_sphere

    ! Local variables

    real (kind=r8) :: redp(1)

    Max_sphere = MAXVAL(fld(:,:,nets:nete))

    redp(1) = Max_sphere
    call pmax_mt(red_max,redp,1,hybrid)
    Max_sphere = red_max%buf(1)

  end function global_maximum

  ! ==========================================================
  ! l1_snorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(8)
  ! for a scalar quantity
  ! ===========================================================

  function l1_snorm(elem,fld,fld_exact,hybrid,npts,nets,nete) result(l1)

    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: fld(npts,npts,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: fld_exact(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l1

    ! Local variables

    real (kind=r8) :: dfld_abs(npts,npts,nets:nete)
    real (kind=r8) :: fld_exact_abs(npts,npts,nets:nete)
    real (kind=r8) :: dfld_abs_int
    real (kind=r8) :: fld_exact_abs_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dfld_abs(i,j,ie) = ABS(fld(i,j,ie)-fld_exact(i,j,ie))
             fld_exact_abs(i,j,ie) = ABS(fld_exact(i,j,ie))
          end do
       end do
    end do

    dfld_abs_int = global_integral(elem, dfld_abs(:,:,nets:nete),hybrid,npts,nets,nete)
    fld_exact_abs_int = global_integral(elem, fld_exact_abs(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dfld_abs_int/fld_exact_abs_int

  end function l1_snorm

  ! ===========================================================
  ! l1_vnorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(97),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l1_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l1)
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l1

    ! Local variables

    real (kind=r8), dimension(:,:,:,:), pointer :: met
    real (kind=r8) :: dvsq(npts,npts,nets:nete)
    real (kind=r8) :: vtsq(npts,npts,nets:nete)
    real (kind=r8) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: dv1,dv2
    real (kind=r8) :: vt1,vt2
    real (kind=r8) :: dvsq_int
    real (kind=r8) :: vtsq_int

    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dvsq_int/vtsq_int

  end function l1_vnorm

  ! ==========================================================
  ! l2_snorm:
  !
  ! computes the l2 norm per Williamson et al, p. 218 eq(83)
  ! for a scalar quantity on the pressure grid.
  !
  ! ===========================================================

  function l2_snorm(elem,fld,fld_exact,hybrid,npts,nets,nete) result(l2)
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t), intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: fld(npts,npts,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: fld_exact(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l2

    ! Local variables

    real (kind=r8) :: dh2(npts,npts,nets:nete)
    real (kind=r8) :: fld_exact2(npts,npts,nets:nete)
    real (kind=r8) :: dh2_int
    real (kind=r8) :: fld_exact2_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dh2(i,j,ie)=(fld(i,j,ie)-fld_exact(i,j,ie))**2
             fld_exact2(i,j,ie)=fld_exact(i,j,ie)**2
          end do
       end do
    end do

    dh2_int = global_integral(elem,dh2(:,:,nets:nete),hybrid,npts,nets,nete)
    fld_exact2_int = global_integral(elem,fld_exact2(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dh2_int)/SQRT(fld_exact2_int)

  end function l2_snorm

  ! ==========================================================
  ! l2_vnorm:
  !
  ! computes the l2 norm per Williamson et al, p. 219 eq(98)
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l2_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l2)
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: l2

    ! Local variables

    real (kind=r8), dimension(:,:,:,:), pointer :: met
    real (kind=r8) :: dvsq(npts,npts,nets:nete)
    real (kind=r8) :: vtsq(npts,npts,nets:nete)
    real (kind=r8) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: dv1,dv2
    real (kind=r8) :: vt1,vt2
    real (kind=r8) :: dvsq_int
    real (kind=r8) :: vtsq_int
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2
             vtsq(i,j,ie) = vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dvsq_int)/SQRT(vtsq_int)

  end function l2_vnorm

  ! ==========================================================
  ! linf_snorm:
  !
  ! computes the l infinity norm per Williamson et al, p. 218 eq(84)
  ! for a scalar quantity on the pressure grid...
  !
  ! ===========================================================

  function linf_snorm(fld,fld_exact,hybrid,npts,nets,nete) result(linf)
    use hybrid_mod, only : hybrid_t
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: fld(npts,npts,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: fld_exact(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: linf

    ! Local variables

    real (kind=r8) :: dfld_abs(npts,npts,nets:nete)
    real (kind=r8) :: fld_exact_abs(npts,npts,nets:nete)
    real (kind=r8) :: dfld_abs_max
    real (kind=r8) :: fld_exact_abs_max
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dfld_abs(i,j,ie)=ABS(fld(i,j,ie)-fld_exact(i,j,ie))
             fld_exact_abs(i,j,ie)=ABS(fld_exact(i,j,ie))
          end do
       end do
    end do

    dfld_abs_max = global_maximum(dfld_abs(:,:,nets:nete),hybrid,npts,nets,nete)
    fld_exact_abs_max = global_maximum(fld_exact_abs(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dfld_abs_max/fld_exact_abs_max

  end function linf_snorm


  ! ==========================================================
  ! linf_vnorm:
  !
  ! computes the linf norm per Williamson et al, p. 218 eq(99),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function linf_vnorm(elem,v,vt,hybrid,npts,nets,nete) result(linf)
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=r8), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=r8), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8)             :: linf

    ! Local variables

    real (kind=r8), dimension(:,:,:,:), pointer :: met
    real (kind=r8) :: dvsq(npts,npts,nets:nete)
    real (kind=r8) :: vtsq(npts,npts,nets:nete)
    real (kind=r8) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=r8) :: dv1,dv2
    real (kind=r8) :: vt1,vt2
    real (kind=r8) :: dvsq_max
    real (kind=r8) :: vtsq_max
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met

       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_max = global_maximum(dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_max = global_maximum(vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dvsq_max/vtsq_max

  end function linf_vnorm

  subroutine wrap_repro_sum (nvars, comm, nsize)
    use dimensions_mod,   only: nelemd
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
    use cam_abortutils,   only: endrun
    use parallel_mod,     only: global_shared_buf, global_shared_sum, nrepro_vars

    integer :: nvars            !  number of variables to be summed (cannot exceed nrepro_vars)
    integer :: comm             !  mpi communicator
    integer, optional :: nsize  !  local buffer size (defaults to nelemd - number of elements in mpi task)

    integer nsize_use

    if (present(nsize)) then
       nsize_use = nsize
    else
       nsize_use = nelemd
    endif
    if (nvars .gt. nrepro_vars) call endrun('ERROR: repro_sum_buffer_size exceeded')

! Repro_sum contains its own OpenMP, so only one thread should call it (AAM)

!$OMP BARRIER
!$OMP MASTER

    call repro_sum(global_shared_buf, global_shared_sum, nsize_use, nelemd, nvars, commid=comm)


!$OMP END MASTER
!$OMP BARRIER

  end subroutine wrap_repro_sum

  subroutine automatically_set_viscosity_coefficients(hybrid,ne,max_min_dx,min_min_dx,nu,factor,str)
    use physconst,      only: rearth
    use control_mod,    only: hypervis_scaling,hypervis_power
    use hybrid_mod,     only: hybrid_t
    use cam_abortutils, only: endrun

    type (hybrid_t), intent(in)    :: hybrid
    integer        , intent(in)    :: ne
    real (kind=r8),  intent(in)    :: max_min_dx,min_min_dx,factor
    real (kind=r8),  intent(inout) :: nu
    character(len=4), intent(in)   :: str

    real(r8)      :: uniform_res_hypervis_scaling,nu_fac
    real(kind=r8) :: nu_min, nu_max
    !
    !************************************************************************************************************
    !
    ! automatically set viscosity coefficients
    !
    !
    ! Use scaling from
    !
    ! - Boville, B. A., 1991: Sensitivity of simulated climate to
    !   model resolution. J. Climate, 4, 469-485.
    !
    ! - TAKAHASHI ET AL., 2006: GLOBAL SIMULATION OF MESOSCALE SPECTRUM
    !
    uniform_res_hypervis_scaling = 1.0_r8/log10(2.0_r8)
    !
    ! compute factor so that at ne30 resolution nu=1E15
    ! scale so that scaling works for other planets
    !
    ! grid spacing in meters = max_min_dx*1000.0_r8
    !
    nu_fac = (rearth/6.37122E6_r8)*1.0E15_r8/(110000.0_r8**uniform_res_hypervis_scaling)

    if (nu < 0) then
      if (ne <= 0) then
        if (hypervis_power/=0) then
          call endrun('ERROR: Automatic scaling of scalar viscosity not implemented')
        else if (hypervis_scaling/=0) then
          nu_min = factor*nu_fac*(max_min_dx*1000.0_r8)**uniform_res_hypervis_scaling
          nu_max = factor*nu_fac*(min_min_dx*1000.0_r8)**uniform_res_hypervis_scaling
          nu     = factor*nu_min
          if (hybrid%masterthread) then
            write(iulog,'(a,a)')             "Automatically setting nu",TRIM(str)
            write(iulog,'(a,2e9.2,a,2f9.2)') "Value at min/max grid spacing: ",nu_min,nu_max,&
                 " Max/min grid spacing (km) = ",max_min_dx,min_min_dx
          end if
          nu = nu_min*(2.0_r8*rearth/(3.0_r8*max_min_dx*1000.0_r8))**hypervis_scaling/(rearth**4)
          if (hybrid%masterthread) &
               write(iulog,'(a,a,a,e9.3)') "Nu_tensor",TRIM(str)," = ",nu
        end if
      else
        nu     = factor*nu_fac*((30.0_r8/ne)*110000.0_r8)**uniform_res_hypervis_scaling
        if (hybrid%masterthread) then
          write(iulog,'(a,a,a,e9.2)') "Automatically setting nu",TRIM(str)," =",nu
        end if
      end if
    end if
  end subroutine automatically_set_viscosity_coefficients
end module global_norms_mod
