!#define fv3remap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Begin GPU remap module  !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vertremap_mod

  !**************************************************************************************
  !
  !  Purpose:
  !        Construct sub-grid-scale polynomials using piecewise spline method with
  !        monotone filters.
  !
  !  References: PCM - Zerroukat et al., Q.J.R. Meteorol. Soc., 2005. (ZWS2005QJR)
  !              PSM - Zerroukat et al., Int. J. Numer. Meth. Fluids, 2005. (ZWS2005IJMF)
  !
  !**************************************************************************************

  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: np,nlev,qsize,nlevp,npsq,nc
  use hybvcoord_mod,          only: hvcoord_t
  use element_mod,            only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use perf_mod,               only: t_startf, t_stopf ! _EXTERNAL
  use parallel_mod,           only: parallel_t
  use cam_abortutils,         only: endrun
  use control_mod,            only: vert_remap_q_alg

  public remap1                  ! remap any field, splines, monotone
  public remap1_nofilter         ! remap any field, splines, no filter
! todo: tweak interface to match remap1 above, rename remap1_ppm:
  public remap_q_ppm             ! remap state%Q, PPM, monotone

  contains

!=======================================================================================================!

subroutine remap1(Qdp,nx,qstart,qstop,qsize,dp1,dp2,identifier,ptop)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
  !
  implicit none
  integer, intent(in) :: nx,qstart,qstop,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  integer, intent(in) :: identifier  !1: tracers, 2: T, 3: u,v
  real (kind=r8), intent(in) :: ptop
  ! ========================
  ! Local Variables
  ! ========================
  real (kind=r8) :: pe1(nx,nlev+1),pe2(nx,nlev+1),inv_dp(nx,nx,nlev),dp2_local(nx,nlev)
  real (kind=r8) :: tmp(nx,nlev), gz(nx)
  integer        :: kord_tr(qsize),i,j,k,itrac
  logical        :: fill
#ifndef fv3remap
  call remap_Q_ppm(qdp,nx,qstart,qstop,qsize,dp1,dp2)
#else
  
  if (identifier==1) then
    !    pe1 interface of Lagrangian levels
    inv_dp = 1.0_r8/dp1
    do itrac=1,qsize
      Qdp(:,:,:,itrac) = Qdp(:,:,:,itrac)*inv_dp(:,:,:)
    end do    
    kord_tr = 9
    !    kord_tr = 4
    fill    = .true.
    do j=1,nx
      pe1(:,1) = ptop
      pe2(:,1) = ptop
      do k=1,nlev
        do i=1,nx
          pe1(i,k+1) = pe1(i,k)+dp1(i,j,k)
          pe2(i,k+1) = pe2(i,k)+dp2(i,j,k)
        end do
      end do
      pe1(:,nlev+1) = pe2(:,nlev+1)
      do k=1,nlev
        do i=1,nx
          dp2_local(i,k) = pe2(i,k+1)-pe2(i,k)
          !          pe2(i,k+1) = pe2(i,k)+dp2(i,j,k)
        end do
      end do
      call mapn_tracer(qsize, nlev, pe1, pe2, Qdp, dp2_local, kord_tr, j,     &
           1, nx, 1, nx, 1, nx, 0.0_r8, fill)
    end do
    
    do itrac=1,qsize
      Qdp(:,:,:,itrac) = Qdp(:,:,:,itrac)*dp2(:,:,:)
    end do
  else if (identifier==2) then
    gz = 0.0_r8 !note: gz not used for iv=-1
    !
    ! mapping of winds
    !
    kord_tr = 9
    do j=1,nx
      pe1(:,1) = ptop
      pe2(:,1) = ptop
      do k=1,nlev
        do i=1,nx
          pe1(i,k+1) = pe1(i,k)+dp1(i,j,k)
          pe2(i,k+1) = pe2(i,k)+dp2(i,j,k)
        end do
      end do
      pe1(:,nlev+1) = pe2(:,nlev+1)

      do k=1,nlev+1
        do i=1,nx
          pe1(i,k) = log(pe1(i,k))
          pe2(i,k) = log(pe2(i,k))
        end do
      end do      
      do itrac=1,qsize
        call map_scalar(nlev,  pe1, Qdp(:,:,:,itrac) , gz,   & 
             nlev,  pe2,          Qdp(:,:,:,itrac) ,              &
             1, nx, j, 1, nx, 1, nx, 1, abs(kord_tr(itrac)), 184.0_r8)        
      end do
    end do
  else if (identifier==3) then
    gz = 0.0_r8 !note: gz not used for iv=-1
    !
    ! mapping of winds
    !
    kord_tr = 9
    do j=1,nx
      pe1(:,1) = ptop
      pe2(:,1) = ptop
      do k=1,nlev
        do i=1,nx
          pe1(i,k+1) = pe1(i,k)+dp1(i,j,k)
          pe2(i,k+1) = pe2(i,k)+dp2(i,j,k)
        end do
      end do
      pe1(:,nlev+1) = pe2(:,nlev+1)
      do itrac=1,qsize
        call map1_ppm( nlev, pe1(:,:),   Qdp(:,:,:,itrac),   gz,   &
             nlev, pe2(:,:),    Qdp(:,:,:,itrac),               &
             1, nx, j, 1, nx, 1, nx, -1, kord_tr(itrac))        
      end do
    end do
  else
    call remap_Q_ppm(qdp,nx,qstart,qstop,qsize,dp1,dp2)
  end if
  
#endif  
end subroutine remap1

subroutine remap1_nofilter(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer, intent(in) :: nx,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=r8), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=r8), dimension(nlev)      :: h,Qcol,za0,za1,za2,zarg,zhdp
  real (kind=r8)  :: tmp_cal,zv1,zv2
  integer :: zkr(nlev+1),i,ilev,j,jk,k,q
  logical :: abort=.false.
  !   call t_startf('remap1_nofilter')

#if (defined COLUMN_OPENMP)
  !$omp parallel do num_threads(tracer_num_threads) &
  !$omp    private(q,i,j,z1c,z2c,zv,k,Qcol,zkr,ilev) &
  !$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal) &
  !$omp    private(za0,za1,za2) &
  !$omp    private(ip2,zv1,zv2)
#endif
  do q=1,qsize
    do i=1,nx
      do j=1,nx

        z1c(1)=0 ! source grid
        z2c(1)=0 ! target grid
        do k=1,nlev
          z1c(k+1)=z1c(k)+dp1(i,j,k)
          z2c(k+1)=z2c(k)+dp2(i,j,k)
        enddo

        zv(1)=0
        do k=1,nlev
          Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
          zv(k+1) = zv(k)+Qcol(k)
        enddo

        if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001_r8) then
          write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
          write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
          write(6,*) 'DATA FOR MODEL LEVELS'
          write(6,*) 'PLEVMODEL=',z2c(nlev+1)
          write(6,*) 'PLEV     =',z1c(nlev+1)
          write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
          abort=.true.
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! quadratic splies with UK met office monotonicity constraints  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        zkr  = 99
        ilev = 2
        zkr(1) = 1
        zkr(nlev+1) = nlev
        kloop: do k = 2,nlev
          do jk = ilev,nlev+1
            if (z1c(jk).ge.z2c(k)) then
              ilev      = jk
              zkr(k)   = jk-1
              cycle kloop
            endif
          enddo
        enddo kloop

        zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
        zgam(1)      = 0.0_r8
        zgam(nlev+1) = 1.0_r8
        zhdp = z1c(2:nlev+1)-z1c(1:nlev)


        h = 1/zhdp
        zarg = Qcol * h
        rhs = 0
        lower_diag = 0
        diag = 0
        upper_diag = 0

        rhs(1)=3*zarg(1)
        rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
        rhs(nlev+1)=3*zarg(nlev)

        lower_diag(1)=1
        lower_diag(2:nlev) = h(1:nlev-1)
        lower_diag(nlev+1)=1

        diag(1)=2
        diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
        diag(nlev+1)=2

        upper_diag(1)=1
        upper_diag(2:nlev) = h(2:nlev)
        upper_diag(nlev+1)=0

        q_diag(1)=-upper_diag(1)/diag(1)
        rhs(1)= rhs(1)/diag(1)

        do k=2,nlev+1
          tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
          q_diag(k) = -upper_diag(k)*tmp_cal
          rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
        enddo
        do k=nlev,1,-1
          rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
        enddo

        za0 = rhs(1:nlev)
        za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
        za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! start iteration from top to bottom of atmosphere !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        zv1 = 0
        do k=1,nlev
          if (zgam(k+1)>1_r8) then
            WRITE(*,*) 'r not in [0:1]', zgam(k+1)
            abort=.true.
          endif
          zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
               (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
          Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
          zv1 = zv2
        enddo
      enddo
    enddo
  enddo ! q loop
  if (abort) then
    call endrun('Bad levels in remap1_nofilter.  usually CFL violatioin')
  end if
end subroutine remap1_nofilter

!=============================================================================!

!#ifndef fv3remap
!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm(Qdp,nx,qstart,qstop,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  use control_mod, only : vert_remap_q_alg

  implicit none
  integer,intent(in) :: nx,qstart,qstop,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=r8), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=r8), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=r8), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=r8), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=r8), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=r8), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=r8), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=r8), dimension(       nlev   ) :: z1, z2
  real(kind=r8) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=r8) :: massn1, massn2, ext(2)
  integer :: i, j, k, q, kk, kid(nlev)

  do j = 1 , nx
    do i = 1 , nx

      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo(k)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo



      pio(nlev+2) = pio(nlev+1) + 1._r8  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo(1   -k) = dpo(       k)
        dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
!        do while ( pio(kk) <= pin(k+1) )
!          kk = kk + 1
!          if(kk==nlev+2) exit
!        enddo
        !        kk = kk - 1                   !kk is now the cell index we're integrating over.

        if (pio(kk) <= pin(k+1)) then
          do while ( pio(kk) <= pin(k+1) )
            kk = kk + 1
          enddo
          kk = kk - 1                   !kk is now the cell index we're integrating over.
        else
          call binary_search(pio, pin(k+1), kk)
        end if
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5_R8                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5_r8 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( dpo )

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = qstart, qstop
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0._r8

        do k = 1 , nlev
          ao(k) = Qdp(i,j,k,q)
          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo
        !Fill in ghost values. Ignored if vert_remap_q_alg == 2
        if (vert_remap_q_alg == 10) then
          ext(1) = minval(ao(1:nlev))
          ext(2) = maxval(ao(1:nlev))
          call linextrap(dpo(2), dpo(1), dpo(0), dpo(-1), ao(2), ao(1), ao(0), ao(-1), ext(1), ext(2))
          call linextrap(dpo(nlev-1), dpo(nlev), dpo(nlev+1), dpo(nlev+2), &
               ao(nlev-1), ao(nlev), ao(nlev+1), ao(nlev+2), ext(1), ext(2))
        else
          do k = 1 , gs
            ao(1   -k) = ao(       k)
            ao(nlev+k) = ao(nlev+1-k)
          enddo
        end if
        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx )
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0._r8
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
          Qdp(i,j,k,q) = massn2 - massn1
          massn1 = massn2
        enddo
      enddo
    enddo
  enddo
! call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm
!#endif
! Find k such that pio(k) <= pivot < pio(k+1). Provide a reasonable input
! value for k.
subroutine binary_search(pio, pivot, k)
  real(kind=r8), intent(in)    :: pio(nlev+2), pivot
  integer,       intent(inout) :: k
  integer :: lo, hi, mid
  
  if (pio(k) > pivot) then
    lo = 1
    hi = k
  else
    lo = k
    hi = nlev+2
  end if
  do while (hi > lo + 1)
    k = (lo + hi)/2
    if (pio(k) > pivot) then
      hi = k
    else
      lo = k
    end if
  end do
  k = lo
end subroutine binary_search
!=======================================================================================================!

!This compute grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids( dx )   result(rslt)
  implicit none
  real(kind=r8), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=r8)             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  do j = 0 , nlev+1
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2._r8*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2._r8*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  do j = 0 , nlev
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1._r8 / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2._r8 * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2._r8 * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2._r8 * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2._r8*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2._r8*dx(j+1) )
  enddo
end function compute_ppm_grids


!=======================================================================================================!



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx )    result(coefs)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=r8), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=r8), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  real(kind=r8) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=r8) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=r8) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real(kind=r8) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=r8) :: al, ar                            !Left and right interface values for cell-local limiting
  integer :: j
  integer :: indB, indE

  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  do j = 0 , nlev+1
    da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
    dma(j) = minval( (/ abs(da) , 2._r8 * abs( a(j) - a(j-1) ) , 2._r8 * abs( a(j+1) - a(j) ) /) ) * sign(1._r8,da)
    if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0._r8 ) dma(j) = 0._r8
  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  do j = 0 , nlev
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
         * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  do j = 1 , nlev
    al = ai(j-1)
    ar = ai(j  )
    if ( (ar - a(j)) * (a(j) - al) <= 0._r8 ) then
      al = a(j)
      ar = a(j)
    endif
    if ( (ar - al) * (a(j) - (al + ar)/2._r8) >  (ar - al)**2/6. ) al = 3._r8*a(j) - 2._r8 * ar
    if ( (ar - al) * (a(j) - (al + ar)/2._r8) < -(ar - al)**2/6. ) ar = 3._r8*a(j) - 2._r8 * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5_r8 * a(j) - ( al + ar ) / 4._r8
    coefs(1,j) = ar - al
    ! coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
    coefs(2,j) = 3._r8 * (-2._r8 * a(j) + ( al + ar ))
  enddo

  !If vert_remap_q_alg == 2, use piecewise constant in the boundaries, and don't use ghost cells.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0._r8
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0._r8
  endif
end function compute_ppm


!=======================================================================================================!


!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=r8), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=r8), intent(in) :: x1      !lower domain bound for integration
  real(kind=r8), intent(in) :: x2      !upper domain bound for integration
  real(kind=r8)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 ** 2 - x1 ** 2) / 0.2D1 + a(2) * (x2 ** 3 - x1 ** 3) / 0.3D1
end function integrate_parabola


!=============================================================================================!
  subroutine linextrap(dx1,dx2,dx3,dx4,y1,y2,y3,y4,lo,hi)
    real(kind=r8), intent(in) :: dx1,dx2,dx3,dx4,y1,y2,lo,hi
    real(kind=r8), intent(out) :: y3,y4

    real(kind=r8), parameter :: half = 0.5_r8

    real(kind=r8) :: x1,x2,x3,x4,a

   x1 = half*dx1
    x2 = x1 + half*(dx1 + dx2)
    x3 = x2 + half*(dx2 + dx3)
    x4 = x3 + half*(dx3 + dx4)
    a  = (x3-x1)/(x2-x1)
    y3 = (1.0_r8-a)*y1 + a*y2
    a  = (x4-x1)/(x2-x1)
    y4 = (1.0_r8-a)*y1 + a*y2
    y3 = max(lo, min(hi, y3))
    y4 = max(lo, min(hi, y4))
  end subroutine linextrap


#ifdef fv3remap
 subroutine map_scalar( km,   pe1,    q1,   qs,           &
                        kn,   pe2,    q2,   i1, i2,       &
                         j,  ibeg, iend, jbeg, jend, iv,  kord, q_min)
! iv=1
 integer, intent(in) :: i1                !< Starting longitude
 integer, intent(in) :: i2                !< Finishing longitude
 integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == temp 2 == remap temp with cs scheme
 integer, intent(in) :: kord              !< Method order
 integer, intent(in) :: j                 !< Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                !< Original vertical dimension
 integer, intent(in) :: kn                !< Target vertical dimension
 real(kind=r8), intent(in) ::   qs(i1:i2)       !< bottom BC
 real(kind=r8), intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
 real(kind=r8), intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
 real(kind=r8), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
! INPUT/OUTPUT PARAMETERS:
 real(kind=r8), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output
 real(kind=r8), intent(in):: q_min

! DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! LOCAL VARIABLES:
   real(kind=r8)    dp1(i1:i2,km)
   real(kind=r8)   q4(4,i1:i2,km)
   real(kind=r8)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map_scalar

  
 subroutine mapn_tracer(nq, km, pe1, pe2, q1, dp2, kord, j,     &
                        i1, i2, isd, ied, jsd, jed, q_min, fill)
! INPUT PARAMETERS:
      integer, intent(in):: km                !< vertical dimension
      integer, intent(in):: j, nq, i1, i2
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: kord(nq)
      real(kind=r8), intent(in)::  pe1(i1:i2,km+1)     !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
      real(kind=r8), intent(in)::  pe2(i1:i2,km+1)     !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
      real(kind=r8), intent(in)::  dp2(i1:i2,km)
      real(kind=r8), intent(in)::  q_min
      logical, intent(in):: fill
      real(kind=r8), intent(inout):: q1(isd:ied,jsd:jed,km,nq) ! Field input
! LOCAL VARIABLES:
      real(kind=r8):: q4(4,i1:i2,km,nq)
      real(kind=r8):: q2(i1:i2,km,nq) !< Field output
      real(kind=r8):: qsum(nq)
      real(kind=r8):: dp1(i1:i2,km)
      real(kind=r8):: qs(i1:i2)
      real(kind=r8):: pl, pr, dp, esl, fac1, fac2
      integer:: i, k, l, m, k0, iq

      do k=1,km
         do i=i1,i2
            dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

      do iq=1,nq
         do k=1,km
            do i=i1,i2
               q4(1,i,k,iq) = q1(i,j,k,iq)
            enddo
         enddo
         call scalar_profile( qs, q4(1,i1,1,iq), dp1, km, i1, i2, 0, kord(iq), q_min )
      enddo

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,km
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            fac1 = pr + pl
            fac2 = r3*(pr*fac1 + pl*pl) 
            fac1 = 0.5*fac1
            do iq=1,nq
               q2(i,k,iq) = q4(2,i,l,iq) + (q4(4,i,l,iq)+q4(3,i,l,iq)-q4(2,i,l,iq))*fac1  &
                                         -  q4(4,i,l,iq)*fac2
            enddo
            k0 = l
            goto 555
          else
! Fractional area...
            dp = pe1(i,l+1) - pe2(i,k)
            fac1 = 1. + pl
            fac2 = r3*(1.+pl*fac1)
            fac1 = 0.5*fac1
            do iq=1,nq
               qsum(iq) = dp*(q4(2,i,l,iq) + (q4(4,i,l,iq)+   &
                              q4(3,i,l,iq) - q4(2,i,l,iq))*fac1 - q4(4,i,l,iq)*fac2)
            enddo
            do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
               if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp1(i,m)*q4(1,i,m,iq)
                  enddo
               else
                  dp = pe2(i,k+1)-pe1(i,m)
                  esl = dp / dp1(i,m)
                  fac1 = 0.5*esl
                  fac2 = 1.-r23*esl
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp*( q4(2,i,m,iq) + fac1*(         &
                                q4(3,i,m,iq)-q4(2,i,m,iq)+q4(4,i,m,iq)*fac2 ) )
                  enddo
                  k0 = m
                  goto 123
               endif
            enddo
            goto 123
          endif
      endif
100   continue
123   continue
      do iq=1,nq
         q2(i,k,iq) = qsum(iq) / dp2(i,k)
      enddo
555   continue
1000  continue

  if (fill) call fillz(i2-i1+1, km, nq, q2, dp2)

  do iq=1,nq
!    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
     do k=1,km
        do i=i1,i2
           q1(i,j,k,iq) = q2(i,k,iq)
        enddo
     enddo
  enddo

 end subroutine mapn_tracer


  subroutine map1_ppm( km,   pe1,    q1,   qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                !< Starting longitude
 integer, intent(in) :: i2                !< Finishing longitude
 integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == ??? 2 == remap temp with cs scheme
 integer, intent(in) :: kord              !< Method order
 integer, intent(in) :: j                 !< Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                !< Original vertical dimension
 integer, intent(in) :: kn                !< Target vertical dimension
 real(kind=r8), intent(in) ::   qs(i1:i2)       !< bottom BC
 real(kind=r8), intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
 real(kind=r8), intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
 real(kind=r8), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
! INPUT/OUTPUT PARAMETERS:
 real(kind=r8), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output

! DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

! LOCAL VARIABLES:
   real(kind=r8)    dp1(i1:i2,km)
   real(kind=r8)   q4(4,i1:i2,km)
   real(kind=r8)    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm

 subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)

! INPUT PARAMETERS:
 integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      !< Starting longitude
 integer, intent(in):: i2      !< Finishing longitude
 integer, intent(in):: km      !< Vertical dimension
 integer, intent(in):: kord    !< Order (or more accurately method no.):
                               ! 
 real(kind=r8) , intent(in):: delp(i1:i2,km)     !< Layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real(kind=r8) , intent(inout):: a4(4,i1:i2,km)  !< Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
! 
! !REVISION HISTORY: 
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real(kind=r8)    dc(i1:i2,km)
      real(kind=r8)    h2(i1:i2,km)
      real(kind=r8)  delq(i1:i2,km)
      real(kind=r8)   df2(i1:i2,km)
      real(kind=r8)    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real(kind=r8)  fac
      real(kind=r8)  a1, a2, c1, c2, c3, d1, d2
      real(kind=r8)  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo 
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

#ifdef BOT_MONO
      do i=i1,i2
         a4(4,i,km) = 0
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
         d1 = a4(1,i,km) - a4(2,i,km)
         d2 = a4(3,i,km) - a4(1,i,km)
         if ( d1*d2 < 0. ) then
              a4(2,i,km) = a4(1,i,km)
              a4(3,i,km) = a4(1,i,km)
         else
              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
              a4(2,i,km) = a4(1,i,km) - dq
              a4(3,i,km) = a4(1,i,km) + dq
         endif
      enddo
#else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
#endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2 
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile

 subroutine ppm_limiters(dm, a4, itot, lmt)

! INPUT PARAMETERS:
      real(kind=r8) , intent(in):: dm(*)     !< Linear slope
      integer, intent(in) :: itot      !< Total Longitudes
      integer, intent(in) :: lmt       !< 0: Standard PPM constraint 1: Improved full monotonicity constraint
                                       !< (Lin) 2: Positive definite constraint 
                                       !< 3: do nothing (return immediately)
! INPUT/OUTPUT PARAMETERS:
      real(kind=r8) , intent(inout) :: a4(4,*)   !< PPM array AA <-- a4(1,i) AL <-- a4(2,i) AR <-- a4(3,i) A6 <-- a4(4,i)
! LOCAL VARIABLES:
      real(kind=r8)  qmp
      real(kind=r8)  da1, da2, a6da
      real(kind=r8)  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters
 
 
 subroutine scalar_profile(qs, a4, delp, km, i1, i2, iv, kord, qmin)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      !< vertical dimension
 integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others
 integer, intent(in):: kord
 real(kind=r8), intent(in)   ::   qs(i1:i2)
 real(kind=r8), intent(in)   :: delp(i1:i2,km)     !< Layer pressure thickness
 real(kind=r8), intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
 real(kind=r8), intent(in):: qmin
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real(kind=r8)  gam(i1:i2,km)
 real(kind=r8)    q(i1:i2,km+1)
 real(kind=r8)   d4(i1:i2)
 real(kind=r8)   bet, a_bot, grat 
 real(kind=r8)   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------
!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints 
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
               q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then 
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. a4(1,i,k)<qmin ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==14 ) then

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==15 ) then   ! Revised abs(kord)=9 scheme
       do i=i1,i2
          if ( ext5(i,k) .and. ext5(i,k-1) ) then
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
          else if ( ext5(i,k) .and. ext5(i,k+1) ) then
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
          else if ( ext5(i,k) .and. a4(1,i,k)<qmin ) then
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
          elseif( ext6(i,k) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11, 13
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1).or.ext5(i,k+1).or.a4(1,i,k)<qmin) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then 
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine scalar_profile

 
 subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      !< vertical dimension
 integer, intent(in):: iv      !< iv =-1: winds
                               !< iv = 0: positive definite scalars
                               !< iv = 1: others
 integer, intent(in):: kord
 real(kind=r8), intent(in)   ::   qs(i1:i2)
 real(kind=r8), intent(in)   :: delp(i1:i2,km)     !< layer pressure thickness
 real(kind=r8), intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real(kind=r8)  gam(i1:i2,km)
 real(kind=r8)    q(i1:i2,km+1)
 real(kind=r8)   d4(i1:i2)
 real(kind=r8)   bet, a_bot, grat 
 real(kind=r8)   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km) 
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo
 
  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif
!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints 
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then 
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==14 ) then

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==15 ) then   ! revised kord=9 scehem
       do i=i1,i2
          if ( ext5(i,k) ) then  ! c90_mp122
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
                  a4(2,i,k) = a4(1,i,k)
                  a4(3,i,k) = a4(1,i,k)
             endif
          elseif( ext6(i,k) ) then
! Check within the smooth region if subgrid profile is non-monotonic
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then 
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile

 subroutine cs_limiters(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real(kind=r8) , intent(inout) :: a4(4,im)   !< PPM array
! LOCAL VARIABLES:
 real(kind=r8)  da1, da2, a6da
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters

!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
 !using PPM instead of splines.
#ifdef xxx 
subroutine remap_Q_ppm(Qdp,nx,qstart,qstop,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  use control_mod, only : vert_remap_q_alg

  implicit none
  integer,intent(in) :: nx,qstart,qstop,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=r8), dimension(np,nlev+1) :: pio    !Pressure at interfaces for old grid
  real(kind=r8), dimension(np,nlev+1) :: pin    !Pressure at interfaces for new grid

  
  real(kind=r8), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=r8), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=r8), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=r8), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=r8), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=r8), dimension(       nlev   ) :: z1, z2
  real(kind=r8) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=r8) :: massn1, massn2, ext(2)
  integer :: i, j, k, q, kk, kid(nlev)

  real(kind=r8) :: gz(np)
  integer, parameter :: kord=9
  !
  ! remap winds
  !
  gz = 0.0_r8
  do j=1,np
      call map1_ppm( nlev, pe0(1:np,:),   u,   gz,   &
                     nlev, pe3(1:np,:),   u,               &
                     1, np, j, 1, np, 1, np , -1, kord_mt)

    end do
    
  end subroutine remap_Q_ppm
#endif

 subroutine fillz(im, km, nq, q, dp)
   integer,  intent(in):: im                !< No. of longitudes
   integer,  intent(in):: km                !< No. of levels
   integer,  intent(in):: nq                !< Total number of tracers
   real(kind=r8) , intent(in)::  dp(im,km)           !< pressure thickness
   real(kind=r8) , intent(inout) :: q(im,km,nq)      !< tracer mixing ratio
! LOCAL VARIABLES:
   logical:: zfix(im)
   real(kind=r8) ::  dm(km)
   integer i, k, ic, k1
   real(kind=r8)  qup, qly, dup, dq, sum0, sum1, fac

   do ic=1,nq
#ifdef DEV_GFS_PHYS
! Bottom up:
      do k=km,2,-1
         k1 = k-1
         do i=1,im
           if( q(i,k,ic) < 0. ) then
               q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
               q(i,k ,ic) = 0.
           endif
         enddo
      enddo
! Top down:
      do k=1,km-1
         k1 = k+1
         do i=1,im
            if( q(i,k,ic) < 0. ) then
                q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
                q(i,k ,ic) = 0.
            endif
         enddo
      enddo
#else
! Top layer
      do i=1,im
         if( q(i,1,ic) < 0. ) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = 0.
          endif
      enddo

! Interior
      zfix(:) = .false.
      do k=2,km-1
         do i=1,im
           if( q(i,k,ic) < 0. ) then
             zfix(i) = .true.
             if ( q(i,k-1,ic) > 0. ) then
! Borrow from above
                dq = min ( q(i,k-1,ic)*dp(i,k-1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k-1,ic) = q(i,k-1,ic) - dq/dp(i,k-1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
             endif
             if ( q(i,k,ic)<0.0 .and. q(i,k+1,ic)>0. ) then
! Borrow from below:
                dq = min ( q(i,k+1,ic)*dp(i,k+1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k+1,ic) = q(i,k+1,ic) - dq/dp(i,k+1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
             endif
          endif
         enddo
      enddo
 
! Bottom layer
      k = km
      do i=1,im
        if( q(i,k,ic)<0. .and. q(i,k-1,ic)>0.) then
             zfix(i) = .true.
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min(qly, qup)
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
             q(i,k,  ic) = q(i,k,  ic) + dup/dp(i,k  )
          endif
      enddo

! Perform final check and non-local fix if needed
      do i=1,im
        if ( zfix(i) ) then
           sum0 = 0.
           do k=2,km
              dm(k) = q(i,k,ic)*dp(i,k)
              sum0 = sum0 + dm(k)
           enddo

           if ( sum0 > 0. ) then
             sum1 = 0.
             do k=2,km
                sum1 = sum1 + max(0., dm(k))
             enddo
             fac = sum0 / sum1
             do k=2,km
                q(i,k,ic) = max(0., fac*dm(k)/dp(i,k))
             enddo
           endif

         endif
      enddo
#endif

   enddo
 end subroutine fillz
  
 
#endif

  
end module vertremap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End GPU remap module    !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
