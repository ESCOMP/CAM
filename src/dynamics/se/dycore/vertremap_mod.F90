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

  implicit none
  
  public remap1                  ! remap any field, splines, monotone
  public remap1_nofilter         ! remap any field, splines, no filter
! todo: tweak interface to match remap1 above, rename remap1_ppm:
  public remap_q_ppm             ! remap state%Q, PPM, monotone

  contains

!=======================================================================================================!

    subroutine remap1(Qdp,nx,qstart,qstop,qsize,dp1,dp2,ptop,identifier,Qdp_mass,kord)
      use fv_mapz,      only: map1_ppm
      ! remap 1 field
      ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
      !         dp1   layer thickness (source)
      !         dp2   layer thickness (target)
      !
      ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
      !
      integer, intent(in) :: nx,qstart,qstop,qsize
      real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
      real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
      integer, intent(in) :: identifier  !0: tracers, 1: T, -1: u,v
      real (kind=r8), intent(in) :: ptop
      logical, intent(in) :: Qdp_mass
      integer, intent(in) :: kord(qsize)
      ! ========================
      ! Local Variables
      ! ========================
      real (kind=r8) :: pe1(nx,nlev+1),pe2(nx,nlev+1),inv_dp(nx,nx,nlev),dp2_local(nx,nlev)
      real (kind=r8) :: tmp(nx,nlev), gz(nx)
      integer        :: i,j,k,itrac
      logical        :: logp
      integer        :: kord_local(qsize)

      kord_local = kord

      if (any(kord(:) >= 0)) then
        if (.not.qdp_mass) then
          do itrac=1,qsize
            if (kord(itrac) >= 0) then            
              Qdp(:,:,:,itrac) = Qdp(:,:,:,itrac)*dp1(:,:,:)
            end if
          end do
        end if        
        call remap_Q_ppm(qdp,nx,qstart,qstop,qsize,dp1,dp2,kord)
        if (.not.qdp_mass) then
          do itrac=1,qsize
            if (kord(itrac) >= 0) then            
              Qdp(:,:,:,itrac) = Qdp(:,:,:,itrac)/dp2(:,:,:)
            end if
          end do
        end if        
      endif
      if (any(kord(:)<0)) then
        !
        ! check if remapping over p or log(p)
        !
        ! can not mix and match here (all kord's must >-20 or <=-20)
        !
        if (any(kord(:)>-20)) then
          kord_local = abs(kord)
          logp    = .false.
        else
          kord_local = abs(kord/10)         
          if (identifier==1) then
            logp    = .true.
          else
            logp    = .false.            
          end if
        end if
        !
        ! modified FV3 vertical remapping
        !        
        if (qdp_mass) then
          inv_dp = 1.0_r8/dp1
          do itrac=1,qsize
            if (kord(itrac)<0) then            
              Qdp(:,:,:,itrac) = Qdp(:,:,:,itrac)*inv_dp(:,:,:)
            end if
          end do
        end if
        if (logp) then
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
              if (kord(itrac)<0) then
                call map1_ppm( nlev, pe1(:,:),   Qdp(:,:,:,itrac),   gz,   &
                     nlev, pe2(:,:),    Qdp(:,:,:,itrac),               &
                     1, nx, j, 1, nx, 1, nx, identifier, kord_local(itrac))
              end if
            end do
            !      call mapn_tracer(qsize, nlev, pe1, pe2, Qdp, dp2_local, kord, j,     &
            !           1, nx, 1, nx, 1, nx, 0.0_r8, fill)
          end do
        else
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
              if (kord(itrac)<0) then
                call map1_ppm( nlev, pe1(:,:),   Qdp(:,:,:,itrac),   gz,   &!phl
                     nlev, pe2(:,:),    Qdp(:,:,:,itrac),               &
                     1, nx, j, 1, nx, 1, nx, identifier, kord_local(itrac))
              end if
            end do
            !      call mapn_tracer(qsize, nlev, pe1, pe2, Qdp, dp2_local, kord, j,     &
            !           1, nx, 1, nx, 1, nx, 0.0_r8, fill)
          end do
        end if
        if (qdp_mass) then
          do itrac=1,qsize
            if (kord(itrac)<0) then
              Qdp(:,:,:,itrac) = Qdp(:,:,:,itrac)*dp2(:,:,:)
            end if
          end do
        end if
      end if
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

        if (ABS(z2c(nlev+1)-z1c(nlev+1)) >= 0.000001_r8) then
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
            if (z1c(jk) >= z2c(k)) then
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

!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm(Qdp,nx,qstart,qstop,qsize,dp1,dp2,kord)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer,intent(in) :: nx,qstart,qstop,qsize
  real (kind=r8), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=r8), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  integer       , intent(in) :: kord(qsize)
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
      ppmdx(:,:) = compute_ppm_grids( dpo)

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = qstart, qstop
        if (kord(q) >= 0) then
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
        !Fill in ghost values. Ignored if kord == 2
        if (kord(q) == 10) then
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
        coefs(:,:) = compute_ppm( ao , ppmdx, kord(q) )
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
        end if
      enddo
    enddo
  enddo
! call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm

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
function compute_ppm( a , dx , kord)    result(coefs)
  implicit none
  real(kind=r8), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=r8), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  integer,       intent(in) :: kord
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
    if ( (ar - al) * (a(j) - (al + ar)/2._r8) >  (ar - al)**2/6._r8 ) al = 3._r8*a(j) - 2._r8 * ar
    if ( (ar - al) * (a(j) - (al + ar)/2._r8) < -(ar - al)**2/6._r8 ) ar = 3._r8*a(j) - 2._r8 * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5_r8 * a(j) - ( al + ar ) / 4._r8
    coefs(1,j) = ar - al
    coefs(2,j) = 3._r8 * (-2._r8 * a(j) + ( al + ar ))
  enddo

  !If kord == 2, use piecewise constant in the boundaries, and don't use ghost cells.
  if (kord == 2) then
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
end module vertremap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End GPU remap module    !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
