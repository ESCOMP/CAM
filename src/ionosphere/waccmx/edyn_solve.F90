module edyn_solve
!
! Prepare stencils and call mudpack PDE solver. This is executed
! by the root task only, following the gather_edyn call in edynamo.F90.
!
  use shr_kind_mod ,only: r8 => shr_kind_r8 ! 8-byte reals
  use cam_logfile  ,only: iulog
  use edyn_params  ,only: finit
  use edyn_maggrid ,only: nmlon,nmlonp1,nmlat,nmlath
  use edyn_maggrid ,only: res_nlev, res_ngrid
  use spmd_utils,   only: masterproc
  use edyn_solver_coefs, only: nc, cee, cofum

  implicit none
  private

  public :: edyn_solve_init
  public :: solve_edyn

!
! Global 2d fields for root task to complete serial part of dynamo.
! The zigmxxx, rhs and rims are gathered from subdomains by in sub
! gather_edyn (edynamo.F90).
!
  real(r8),allocatable, dimension(:,:), public :: &
    zigm11_glb  ,&
    zigm22_glb  ,&
    zigmc_glb   ,&
    zigm2_glb   ,&
    rhs_glb
  real(r8),allocatable, dimension(:,:,:), public :: &
    rim_glb    ! pde solver output
  real(r8),allocatable, dimension(:,:) :: &
    phisolv
!
! Dimensions of the grid resolutions for the multi-grid PDE:
  integer, public, protected :: &
    nmlon0, &
    nmlat0, &
    nmlon1, &
    nmlat1, &
    nmlon2, &
    nmlat2, &
    nmlon3, &
    nmlat3, &
    nmlon4, &
    nmlat4, &
    nmlon5, &
    nmlat5, &
    nmlon6, &
    nmlat6, &
    nmlon7, &
    nmlat7
!
! Space needed for descretized coefficients of of dynamo pde at all levels:
!
  integer :: ncee
!
! The following parameters nc0,nc1,... are pointers to the beginning of
!   the coefficients for each level of resolution.
!
  integer :: &
    nc0, &
    nc1, &
    nc2, &
    nc3, &
    nc4, &
    nc5, &
    nc6, &
    nc7

  real(r8), private, pointer :: &
    c0(:),  &
    c1(:),  &
    c2(:),  &
    c3(:),  &
    c4(:),  &
    c5(:),  &
    c6(:),  &
    c7(:)

! phihm is high-latitude potential, set by the high-latitude potential model (e.g. Heelis)
! or is prescribed (e.g. AMIE)
!
  real(r8), allocatable, public :: phihm(:,:) ! high-latitude potential
  real(r8), allocatable, public :: pfrac(:,:) ! NH fraction of potential

  contains

!-----------------------------------------------------------------------
  subroutine edyn_solve_init
    use infnan, only: nan, assignment(=)

    allocate(zigm11_glb(nmlonp1,nmlat))
    allocate(zigm22_glb(nmlonp1,nmlat))
    allocate(zigmc_glb(nmlonp1,nmlat))
    allocate(zigm2_glb(nmlonp1,nmlat))
    allocate(rhs_glb(nmlonp1,nmlat))
    allocate(rim_glb(nmlonp1,nmlat,2))
    allocate(phisolv(0:nmlonp1,0:nmlat+1))

    nmlon0=nmlon+1
    nmlat0=(nmlat +1)/2
    nmlon1=(nmlon0+1)/2
    nmlat1=(nmlat0+1)/2
    nmlon2=(nmlon1+1)/2
    nmlat2=(nmlat1+1)/2
    nmlon3=(nmlon2+1)/2
    nmlat3=(nmlat2+1)/2
    nmlon4=(nmlon3+1)/2
    nmlat4=(nmlat3+1)/2
    nmlon5=(nmlon4+1)/2
    nmlat5=(nmlat4+1)/2
    nmlon6=(nmlon5+1)/2
    nmlat6=(nmlat5+1)/2
    nmlon7=(nmlon6+1)/2
    nmlat7=(nmlat6+1)/2

    allocate(cofum(nmlon0,nmlat0,9))

    ncee=10*nmlon0*nmlat0+9*(nmlon1*nmlat1+nmlon2*nmlat2+nmlon3* &
         nmlat3+nmlon4*nmlat4+nmlon5*nmlat5+nmlon6*nmlat6+nmlon7*nmlat7)

    allocate(cee(ncee))

    nc0=1
    nc1=nc0+10*nmlon0*nmlat0
    nc2=nc1+9 *nmlon1*nmlat1
    nc3=nc2+9 *nmlon2*nmlat2
    nc4=nc3+9 *nmlon3*nmlat3
    nc5=nc4+9 *nmlon4*nmlat4
    nc6=nc5+9 *nmlon5*nmlat5
    nc7=nc6+9 *nmlon6*nmlat6

    c0 => cee
    c1 => cee(nc1:)
    c2 => cee(nc2:)
    c3 => cee(nc3:)
    c4 => cee(nc4:)
    c5 => cee(nc5:)
    c6 => cee(nc6:)
    c7 => cee(nc7:)

    allocate(phihm(nmlonp1,nmlat))
    allocate(pfrac(nmlonp1,nmlat0))

    phihm = nan
    pfrac = nan

  end subroutine edyn_solve_init

!-----------------------------------------------------------------------
  subroutine solve_edyn
!
! Set up stencils for solver:
!
    call stencils
!
! Call mudpack PDE solver:
!
    call solver(cofum,c0)

  end subroutine solve_edyn
!-----------------------------------------------------------------------
  subroutine stencils
  use edyn_params ,only: pi_dyn
  use edyn_maggrid,only: dlatm,dlonm
!
! Locals:
    integer :: i,j,jj,jjj,j0,n,ncc,nmaglon,nmaglat, ndx1,ndx2
    real(r8) :: sym
    real(r8) :: cs(nmlat0)

!
! Set index array nc and magnetic latitude cosine array:
! nc pointes to the start of the coefficient array for each level
    nc(1) = nc0
    nc(2) = nc1
    nc(3) = nc2
    nc(4) = nc3
    nc(5) = nc4
    nc(6) = nc5
    nc(7) = nc6
    nc(8) = nc7
    nc(9) = ncee

    do j=1,nmlat0
      cs(j) = cos(pi_dyn/2._r8-(nmlat0-j)*dlatm)
    enddo ! j=1,nmlat0
!
! Set up difference coefficients. Replace zigm11 by A, zigm22 by B,
! zigmc by C, and zigm2 by D.
!
    j0 = nmlat0-nmlath
    do j=1,nmlath       !  1,49 (assuming nmlat=97)
      jj = nmlath+j-1   ! 49,97
      jjj = nmlath-j+1  ! 49,1
!
! factor 4 from 5-point diff. stencil
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
! Sigma_(lam lam)*cos(lam_m)*DT0DTS/(Delta lam)^2
! -zigmc_north = southern hemis. 49,1 equator-pole
! -zigm2_north = southern hemis. 49,1 equator-pole
!  zigm22 = southern hemis. 49,1 equator-pole
!
      do i=1,nmlonp1
        zigmc_glb(i,jj) = (zigmc_glb(i,jj)+zigm2_glb(i,jj))/ &
          (4._r8*dlatm*dlonm)
        zigm2_glb(i,jj) = zigmc_glb(i,jj)-2._r8*zigm2_glb(i,jj)/ &
          (4._r8*dlatm*dlonm)
        zigm22_glb(i,jj)  = zigm22_glb(i,jj)*cs(j0+j)/dlatm**2
        zigmc_glb(i,jjj)  = -zigmc_glb(i,jj)
        zigm2_glb(i,jjj)  = -zigm2_glb(i,jj)
        zigm22_glb(i,jjj) = zigm22_glb(i,jj)
      enddo ! i=1,nmlonp1
      if (j /= nmlath) then
!
! Sigma_(phi phi)/( cos(lam_m)*DT0DTS*(Delta lon)^2 )
! zigm11 = southern hemis. 49,1 equator-pole
!
        do i = 1,nmlonp1
          zigm11_glb(i,jj)  = zigm11_glb(i,jj)/(cs(j0+j)*dlonm**2)
          zigm11_glb(i,jjj) = zigm11_glb(i,jj)
        enddo
      endif
    enddo ! j=1,nmlath
!
! Set zigm11 to zero at megnetic poles to avoid floating exception
! (values at poles are not used):
!
    do i = 1,nmlonp1
      zigm11_glb(i,1) = 0.0_r8
      zigm11_glb(i,nmlat) = 0.0_r8
    enddo

! Clear array for difference stencils at all levels:
    call clearcee(cee,nmlon0,nmlat0)
!
! Init cofum coefficients:
    cofum(:,:,:) = finit
!
! Calculate contribution to stencils from each PDE coefficient
!
! Sigma_(phi phi)/( cos(lam_m)*dt0dts*(Delta lon)^2 )
    sym = 1._r8
    call stencmd(zigm11_glb,nmlon0,nmlat0,sym,cee,1)
!
! Sigma_(lam lam)*cos(lam_m)*dt0dts/(Delta lam)^2
    sym = 1._r8
    call stencmd(zigm22_glb,nmlon0,nmlat0,sym,cee,4)
!
! Sigma_(phi lam)/( 4*Delta lam* Delta lon )
    sym = -1._r8
    call stencmd(zigmc_glb,nmlon0,nmlat0,sym,cee,2)
!
! Sigma_(lam phi)/( 4*Delta lam* Delta lon )
    sym = -1._r8
    call stencmd(zigm2_glb,nmlon0,nmlat0,sym,cee,3)
!
! Insert RHS in finest stencil:
    do j = 1,nmlat0
      jj = nmlath-nmlat0+j
      do i = 1,nmlon0
        ndx1 = 9*nmlat0*nmlon0 + (j-1)*nmlon0 + i
        c0(ndx1) = rhs_glb(i,jj)
      enddo ! i = 1,nmlon0
    enddo ! j = 1,nmlat0
    ndx1 = 9*nmlat0*nmlon0 + nmlonp1
    ndx2 = 9*nmlat0*nmlon0 + 1
    c0(ndx1) = c0(ndx2)
!
! Set boundary condition at the pole:
    call edges(c0,nmlon0,nmlat0)
    call edges(c1,nmlon1,nmlat1)
    call edges(c2,nmlon2,nmlat2)
    call edges(c3,nmlon3,nmlat3)
    call edges(c4,nmlon4,nmlat4)
    if ( res_nlev > 5 ) then
       call edges(c5,nmlon5,nmlat5)
    endif
    if ( res_nlev > 6 ) then
       call edges(c6,nmlon6,nmlat6)
    endif
    if ( res_nlev > 7 ) then
       call edges(c7,nmlon7,nmlat7)
    endif
    call edges(cofum,nmlon0,nmlat0)
!
! Divide stencils by cos(lam_0) (not rhs):
    call divide(c0,nmlon0,nmlat0,nmlon0,cs,1)
    call divide(c1,nmlon1,nmlat1,nmlon0,cs,1)
    call divide(c2,nmlon2,nmlat2,nmlon0,cs,1)
    call divide(c3,nmlon3,nmlat3,nmlon0,cs,1)
    call divide(c4,nmlon4,nmlat4,nmlon0,cs,1)
    if ( res_nlev > 5 ) then
       call divide(c5,nmlon5,nmlat5,nmlon0,cs,1)
    endif
    if ( res_nlev > 6 ) then
       call divide(c6,nmlon6,nmlat6,nmlon0,cs,1)
    endif
    if ( res_nlev > 7 ) then
       call divide(c7,nmlon7,nmlat7,nmlon0,cs,1)
    endif
    call divide(cofum,nmlon0,nmlat0,nmlon0,cs,0)
!
! Set value of solution to 1. at pole:
    do i=1,nmlon0
      ndx1 = 9*nmlat0*nmlon0 + (nmlat0-1)*nmlon0 + i
      c0(ndx1) = 1._r8
    enddo
!
! Modify stencils and RHS so that the NH high lat potential is inserted at
!  high latitude.  The SH high lat potential will be added back later.
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
! The dynamo is symmetric about the magnetic equator, but the high latitude
!  is anti-symmetric in both hemispheres.  However, since Mudpack uses the
!  NH potential pattern, then the SH potential pattern must be added
!  back into the 2-D phim before the call threed, and before it is
!  transformed to geographic coordinates.
!
    ncc = 1
    nmaglon = nmlon0
    nmaglat = nmlat0
    do n=1,res_nlev ! resolution levels
      call stenmd(nmaglon,nmaglat,cee(ncc),phihm(1,nmlat0),pfrac)
      ncc = ncc+9*nmaglon*nmaglat
      if (n==1) ncc = ncc+nmaglon*nmaglat ! rhs is in 10th slot
      nmaglon = (nmaglon+1)/2
      nmaglat = (nmaglat+1)/2
    enddo

  end subroutine stencils
!-----------------------------------------------------------------------
  subroutine clearcee(cee,nlon0,nlat0)
!
! Zero C arrays for stencil coefficients.
! Cee will contain:
!   c0(nmlon0,nmlat0,10), c1(nmlon1,nmlat1,9), c2(nmlon2,nmlat2,9),
!   c3(nmlon3,nmlat3,9),  c4(nmlon4,nmlat4,9)
!
! Args:
    integer,intent(in) :: nlon0,nlat0
    real(r8),intent(out) :: cee(*)
!
! Local:
    integer :: nlon,nlat,n,m,i
!
! Compute total size of cee
    nlon = nlon0
    nlat = nlat0
    n = 0
    do m=1,res_nlev ! resolution levels
      n = n+nlon*nlat
      nlon = (nlon+1)/2
      nlat = (nlat+1)/2
    enddo
    n = 9*n+nlon0*nlat0
!
! Clear cee:
    do i=1,n
      cee(i) = 0._r8
    enddo
  end subroutine clearcee
!-----------------------------------------------------------------------
  subroutine stencmd(zigm,nlon0,nlat0,sym,cee,ncoef)
!
! Calculate contribution fo 3 by 3 stencil from coefficient zigm
! at each grid point and level.
!
! Args:
    integer,intent(in) ::  &
      nlon0,               & ! longitude dimension of finest grid level
      nlat0,               & ! latitude dimension of finest grid level
      ncoef                  ! integer identifier of coefficient
    real(r8),intent(in) :: &
      zigm(nlon0,nlat0),   & ! coefficients (nlon0+1/2,(nlat0+1)/2)
      sym                    !  1. if zigm symmetric w.r.t. equator, -1 otherwise
    real(r8),intent(inout) ::  & ! output stencil array consisting of c0,c1,c2,c3,c4
      cee(*)
!
! Local:
    integer :: nc,nlon,nlat,n
    real(r8) :: wkarray(-res_ngrid+1:nmlon0+res_ngrid,nmlat0)
!
! Perform half-way interpolation and extend zigm in wkarray:
!
    call htrpex(zigm,nlon0,nlat0,sym,wkarray)
!
! Calculate contribution to stencil for each grid point and level:
!
    nc = 1
    nlon = nlon0
    nlat = nlat0
!
! Calculate modified and unmodified stencil on finest grid
!
    call cnmmod(nlon0,nlon,nlat,cee(nc),ncoef,wkarray,cofum)
!
! Stencils on other grid levels remain the same.
    nc = nc+10*nlon*nlat
    nlon = (nlon+1)/2
    nlat = (nlat+1)/2
!
    do n=2,res_nlev
      call cnm(nlon0,nlon,nlat,cee(nc),ncoef,wkarray)
      nc = nc+9*nlon*nlat
      if (n==1) nc = nc+nlon*nlat
      nlon = (nlon+1)/2
      nlat = (nlat+1)/2
    enddo ! n=1,5
  end subroutine stencmd
!-----------------------------------------------------------------------
  subroutine htrpex(coeff,nmlon0,nmlat0,sym,wkarray)
!
! Perform half-way interpolation on array coeff and extend over 16 grid
! points. Result returned in wkarray.
!
! Args:
    integer,intent(in)   :: nmlon0,nmlat0
    real(r8),intent(in)  :: coeff(nmlon0,nmlat0),sym
    real(r8),intent(out) :: wkarray(-res_ngrid+1:nmlon0+res_ngrid,nmlat0)
!
! Local:
    integer :: i,j,jj
!
! Copy coeff into positions in wkarray:
    do j=1,nmlat0
      jj = nmlat0-j+1
      do i=1,nmlon0
        wkarray(i,j) = sym*coeff(i,jj)
      enddo ! i=1,nmlon0
    enddo ! j=1,nmlat0
!
! Extend over 2*res_ngrid grid spaces to allow for a total of res_nlev grid levels:
    do i=1,res_ngrid
      do j=1,nmlat0
        wkarray(1-i,j) = wkarray(nmlon0-i,j)
        wkarray(nmlon0+i,j) = wkarray(1+i,j)
      enddo ! j=1,nmlat0
    enddo ! i=1,res_ngrid
  end subroutine htrpex
!-----------------------------------------------------------------------
  subroutine cnm(nlon0,nlon,nlat,c,ncoef,wkarray)
!
! Compute contribution to stencil from zigm(ncoef) on grid nlon by nlat,
! Finest grid is nlon0.
!
! Args:
    integer,intent(in) :: &
      nlon0,              & ! finest grid dimensions
      nlon,nlat             ! output grid dimensions
    real(r8),intent(in) :: wkarray(-res_ngrid+1:nmlon0+res_ngrid,nmlat0)
!
! ncoef: integer id of coefficient:
! ncoef = 1 for zigm11
! ncoef = 2 for zigm12 (=zigmc+zigm2)
! ncoef = 3 for zigm21 (=zigmc-zigm2)
! ncoef = 4 for zigm22
!
    integer,intent(in) :: ncoef
    real(r8),intent(inout) :: &
      c(nlon,nlat,*) ! output array for grid point stencils at resolution nlon x nlat
!
! Local:
    integer :: i,j,nint,i0,j0
! For now, retain this pi to insure bit compatability w/ old code
    real(r8),parameter :: pi=3.141592654_r8
    real(r8) :: wk(nlon0,3)
!
! Compute separation of grid points of resolution nlon x nlat within
! grid of resolution nlon0. Evaluate dlon and dlat, grid spacing
! of nlon x nlat.
!
    nint = (nlon0-1)/(nlon-1)
!
! Scan wkarray nlon x nlat calculating and adding contributions to stencil
! from zigm(ncoef)
    i0 = 1-nint
    j0 = 1-nint
!
! zigm11:
! am 2001-6-27 include boundary condition at equator
    if (ncoef==1) then
      do j = 1,nlat-1
        do i = 1,nlon
          c(i,j,1) = c(i,j,1)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i+1)*nint,j0+j*nint))
          c(i,j,5) = c(i,j,5)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          c(i,j,9) = c(i,j,9)-0.5_r8*(wkarray(i0+(i+1)*nint,j0+j*nint)+ &
            2._r8*wkarray(i0+i*nint,j0+j*nint)+wkarray(i0+(i-1)*nint,j0+j*nint))
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! zigm12 (=zigmc+zigm2)
    elseif (ncoef==2) then
      do j = 2,nlat-1
        do i = 1,nlon
          c(i,j,2) = c(i,j,2)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i+1)*nint,j0+j*nint))
          c(i,j,4) = c(i,j,4)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          c(i,j,6) = c(i,j,6)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          c(i,j,8) = c(i,j,8)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i+1)*nint,j0+j*nint))
          wk(i,1) = 0.5_r8*(wkarray(i0+(i+1)*nint,j0+j*nint)- &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          wk(i,2) = (c(i,j,3)+wk(i,1))*(c(i,j,7)-wk(i,1))
          wk(i,3) = sign(wk(i,1),c(i,j,3)+c(i,j,7))
          if (wk(i,2) >= 0._r8) wk(i,3) = 0._r8
          c(i,j,3) = c(i,j,3)+wk(i,1)+wk(i,3)
          c(i,j,7) = c(i,j,7)-wk(i,1)+wk(i,3)
          c(i,j,9) = c(i,j,9)-2._r8*wk(i,3)
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! zigm21 (=zigmc-zigm2)
    elseif (ncoef==3) then
      do j = 2,nlat-1
        do i = 1,nlon
          c(i,j,2) = c(i,j,2)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,4) = c(i,j,4)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,6) = c(i,j,6)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          c(i,j,8) = c(i,j,8)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          wk(i,1) = 0.5_r8*(wkarray(i0+i*nint,j0+(j+1)*nint)- &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
          wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
          if (wk(i,2) >= 0._r8) wk(i,3) = 0._r8
          c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
          c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
          c(i,j,9) = c(i,j,9)-2._r8*wk(i,3)
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
      j = 1
      do i=1,nlon
        c(i,j,2) = c(i,j,2)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))
        c(i,j,4) = c(i,j,4)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))
        wk(i,1) = 0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))
        wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
        wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
        if (wk(i,2) >= 0._r8) wk(i,3) = 0._r8
        c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
        c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
        c(i,j,9) = c(i,j,9)-2._r8*wk(i,3)
      enddo ! i=1,nlon
!
! zigm22:
    elseif (ncoef==4) then
      do j = 2,nlat-1
        do i = 1,nlon
          c(i,j,3) = c(i,j,3)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,7) = c(i,j,7)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          c(i,j,9) = c(i,j,9)-0.5_r8*(wkarray(i0+i*nint,j0+(j-1)*nint)+ &
            2._r8*wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
      j = 1
      do i=1,nlon
        c(i,j,3) = c(i,j,3)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))
        c(i,j,9) = c(i,j,9)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
                                wkarray(i0+i*nint,j0+(j+1)*nint))
      enddo ! i=1,nlon
    endif ! ncoef
  end subroutine cnm
!-----------------------------------------------------------------------
  subroutine cnmmod(nlon0,nlon,nlat,c,ncoef,wkarray,cofum)
!
! Compute contribution to stencil from zigm(ncoef) on grid nlon by nlat,
! Finest grid is nlon0.
!
! Args:
    integer,intent(in) :: &
      nlon0,              & ! finest grid dimensions
      nlon,nlat             ! output grid dimensions
    real(r8),intent(in) :: wkarray(-res_ngrid+1:nmlon0+res_ngrid,nmlat0)
    real(r8),dimension(nmlon0,nmlat0,9),intent(inout) :: cofum
!
! ncoef: integer id of coefficient:
! ncoef = 1 for zigm11
! ncoef = 2 for zigm12 (=zigmc+zigm2)
! ncoef = 3 for zigm21 (=zigmc-zigm2)
! ncoef = 4 for zigm22
!
    integer,intent(in) :: ncoef
    real(r8),intent(inout) :: &
      c(nlon,nlat,*)  ! output array for grid point stencils at resolution nlon x nlat
!
! Local:
    integer :: i,j,nint,i0,j0
! For now, retain this pi to insure bit compatability w/ old code
    real(r8),parameter :: pi=3.141592654_r8
    real(r8) :: wk(nlon0,3)
!
! Compute separation of grid points of resolution nlon x nlat within
! grid of resolution nlon0. Evaluate dlon and dlat, grid spacing
! of nlon x nlat.
!
    nint = (nlon0-1)/(nlon-1)
!
! Scan wkarray nlon x nlat calculating and adding contributions to stencil
! from zigm(ncoef)
    i0 = 1-nint
    j0 = 1-nint
!
! zigm11:
! am 2001-6-27 include boundary condition at equator
    if (ncoef==1) then
      do j = 1,nlat-1
        do i = 1,nlon
          c(i,j,1) = c(i,j,1)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i+1)*nint,j0+j*nint))
          c(i,j,5) = c(i,j,5)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          c(i,j,9) = c(i,j,9)-0.5_r8*(wkarray(i0+(i+1)*nint,j0+j*nint)+ &
            2._r8*wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
!
! Unmodified:
          cofum(i,j,1) = cofum(i,j,1)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i+1)*nint,j0+j*nint))
          cofum(i,j,5) = cofum(i,j,5)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          cofum(i,j,9) = cofum(i,j,9)-0.5_r8*(wkarray(i0+(i+1)*nint,j0+j*nint)+ &
            2._r8*wkarray(i0+i*nint,j0+j*nint)+wkarray(i0+(i-1)*nint,j0+j*nint))
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! zigm12 (=zigmc+zigm2)
    elseif (ncoef==2) then
      do j = 2,nlat-1
        do i = 1,nlon
          c(i,j,2) = c(i,j,2)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i+1)*nint,j0+j*nint))
          c(i,j,4) = c(i,j,4)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          c(i,j,6) = c(i,j,6)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i-1)*nint,j0+j*nint))
          c(i,j,8) = c(i,j,8)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+(i+1)*nint,j0+j*nint))
          wk(i,1) = 0.5_r8*(wkarray(i0+(i+1)*nint,j0+j*nint)- &
            wkarray(i0+(i-1)*nint,j0+j*nint))
!
! Unmodified:
          cofum(i,j,2) = c(i,j,2)
          cofum(i,j,4) = c(i,j,4)
          cofum(i,j,6) = c(i,j,6)
          cofum(i,j,8) = c(i,j,8)
          cofum(i,j,3) = cofum(i,j,3)+wk(i,1)
          cofum(i,j,7) = cofum(i,j,7)-wk(i,1)
!
          wk(i,2) = (c(i,j,3)+wk(i,1))*(c(i,j,7)-wk(i,1))
          wk(i,3) = sign(wk(i,1),c(i,j,3)+c(i,j,7))
          if (wk(i,2) >= 0._r8) wk(i,3) = 0._r8
          c(i,j,3) = c(i,j,3)+wk(i,1)+wk(i,3)
          c(i,j,7) = c(i,j,7)-wk(i,1)+wk(i,3)
          c(i,j,9) = c(i,j,9)-2._r8*wk(i,3)
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! zigm21 (=zigmc-zigm2)
    elseif (ncoef==3) then
      do j = 2,nlat-1
        do i = 1,nlon
          c(i,j,2) = c(i,j,2)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,4) = c(i,j,4)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,6) = c(i,j,6)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          c(i,j,8) = c(i,j,8)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          wk(i,1) = 0.5_r8*(wkarray(i0+i*nint,j0+(j+1)*nint)- &
            wkarray(i0+i*nint,j0+(j-1)*nint))
!
! Unmodified:
          cofum(i,j,2) = c(i,j,2)
          cofum(i,j,4) = c(i,j,4)
          cofum(i,j,6) = c(i,j,6)
          cofum(i,j,8) = c(i,j,8)
          cofum(i,j,1) = cofum(i,j,1)+wk(i,1)
          cofum(i,j,5) = cofum(i,j,5)-wk(i,1)
!
          wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
          wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
          if (wk(i,2) >= 0._r8) wk(i,3) = 0._r8
          c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
          c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
          c(i,j,9) = c(i,j,9)-2._r8*wk(i,3)
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
      j = 1
      do i=1,nlon
        c(i,j,2) = c(i,j,2)+.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))
        c(i,j,4) = c(i,j,4)-.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))
        wk(i,1) = .5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))

        cofum(i,j,2) = c(i,j,2)
        cofum(i,j,4) = c(i,j,4)
        cofum(i,j,1) = cofum(i,j,1)+wk(i,1)
        cofum(i,j,5) = cofum(i,j,5)-wk(i,1)

        wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
        wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
        if (wk(i,2) >= 0._r8) wk(i,3) = 0._r8
        c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
        c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
        c(i,j,9) = c(i,j,9)-2._r8*wk(i,3)
      enddo ! i=1,nlon
!
! zigm22:
    elseif (ncoef==4) then
      do j = 2,nlat-1
        do i = 1,nlon
          c(i,j,3) = c(i,j,3)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,7) = c(i,j,7)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          c(i,j,9) = c(i,j,9)-0.5_r8*(wkarray(i0+i*nint,j0+(j-1)*nint)+ &
            2._r8*wkarray(i0+i*nint,j0+j*nint)+wkarray(i0+i*nint,j0+(j+1)*nint))
!
! Unmodified:
          cofum(i,j,3) = cofum(i,j,3)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,7) = cofum(i,j,7)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
            wkarray(i0+i*nint,j0+(j-1)*nint))
          cofum(i,j,9) = cofum(i,j,9)-0.5_r8*(wkarray(i0+i*nint,j0+(j-1)*nint)+ &
            2._r8*wkarray(i0+i*nint,j0+j*nint)+wkarray(i0+i*nint,j0+(j+1)*nint))
        enddo ! i = 1,nlon
      enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
      j = 1
      do i=1,nlon
        c(i,j,3) = c(i,j,3)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
          wkarray(i0+i*nint,j0+(j+1)*nint))
        c(i,j,9) = c(i,j,9)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
                                wkarray(i0+i*nint,j0+(j+1)*nint))
        cofum(i,j,3) = cofum(i,j,3)+0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
                                wkarray(i0+i*nint,j0+(j+1)*nint))
        cofum(i,j,9) = cofum(i,j,9)-0.5_r8*(wkarray(i0+i*nint,j0+j*nint)+ &
                                wkarray(i0+i*nint,j0+(j+1)*nint))
      enddo ! i=1,nlon
    endif ! ncoef
  end subroutine cnmmod
!--------------------------------------------------------------------
  subroutine edges(c,nlon,nlat)
!
! Insert equatorial and polar boundary conditions in stencil c(nlon,nlat,9)
!
! Args:
    integer,intent(in) :: nlon,nlat
    real(r8),intent(out) :: c(*)
!
! Local:
    integer :: n,i, ndx

    do n=1,8
       do i=1,nlon
        ndx = (n-1)*nlat*nlon + (nlat-1)*nlon + i
        c(ndx) = 0._r8
      enddo
    enddo
    do i=1,nlon
      ndx = 8*nlat*nlon + (nlat-1)*nlon + i
      c(ndx) = 1._r8
    enddo
  end subroutine edges
!--------------------------------------------------------------------
  subroutine divide(c,nlon,nlat,nlon0,cs,igrid)
!
! Divide stencil C by cos(theta(i,j))
!
! Args:
    integer,intent(in) :: nlon,nlat,nlon0,igrid
    real(r8),intent(in) :: cs(:)
    real(r8),intent(out) :: c(*)
!
! Local:
    integer :: nint,j0,n,j,i, ndx
!
    nint = (nlon0-1)/(nlon-1)
    j0 = 1-nint
    do n = 1,9
      do j = 1,nlat-1
        do i = 1,nlon
          ndx = (n-1)*nlat*nlon + (j-1)*nlon + i
          c(ndx) = c(ndx)/(cs(j0+j*nint)*nint**2)
        enddo ! i = 1,nlon
      enddo ! j = 1,nlat-1
    enddo ! n = 1,9
!
    if (nint==1.and.igrid > 0) then
      do i = 1,nlon
        ndx = 9*nlat*nlon + i
        c(ndx) = c(ndx)/cs(1)
      enddo ! i = 1,nlon
    endif
  end subroutine divide
!--------------------------------------------------------------------
  subroutine stenmd(inlon,inlat,c,phihm,pfrac)
  use edyn_params ,only: dtr
  use edyn_maggrid,only: dlatm
!
! Modify stencil to set potential to heelis value within auroral circle.
!
! Args:
    integer,intent(in) :: inlon,inlat
    real(r8),intent(inout) :: c(inlon,inlat,*)
    real(r8),dimension(nmlon0,nmlat0),intent(in) :: &
      phihm,  & ! heelis potential (from subs potm, flwv32)
      pfrac     ! fractional presence of dynamo (from sub colath)
!
! Local:
    integer :: nint,i0,j0,i,j,n,jj
    real(r8) :: real8
!
! Compute separation of grid points for this resolution:
    nint = (nmlon0-1)/(inlon-1)
    i0 = 1-nint
    j0 = 1-nint
!
! If nint==1, then we are at the highest resolution.
! Correct RHS, which is in c(10)
!
    if (nint==1) then
      do j=1,inlat
        do i=1,inlon
          c(i,j,10) = pfrac(i,j)*c(i,j,10)+(1._r8-pfrac(i,j))*c(i,j,9)* &
            (dlatm/(10._r8*dtr))**2*phihm(i,j)
        enddo ! i=1,inlon
      enddo ! j=1,inlat
    endif
!
! Modify stencil, c(i,j,n),n=1,9:
!
    real8 = dble(nint)
    if (nint==1) then
      do j=1,inlat
        jj = j0+j*nint
        do n = 1,8
          do i = 1,inlon
            c(i,j,n) = c(i,j,n)*pfrac(i0+i*nint,jj)
            cofum(i,j,n) = cofum(i,j,n)*pfrac(i0+i*nint,jj)
          enddo ! i = 1,inlon
        enddo ! n = 1,8
        do i = 1,inlon
          c(i,j,9) = c(i,j,9)*pfrac(i0+i*nint,jj)+ &
            (1._r8-pfrac(i0+i*nint,jj))*c(i,j,9)* &
            (dlatm*real8/(10._r8*dtr))**2
          cofum(i,j,9) =cofum(i,j,9)*pfrac(i0+i*nint,jj)+ &
            (1._r8-pfrac(i0+i*nint,jj))*cofum(i,j,9)* &
            (dlatm*real8/(10._r8*dtr))**2
        enddo ! i = 1,inlon
      enddo ! j=1,inlat
    else ! nint /= 1
      do j=1,inlat
        jj = j0+j*nint
        do n = 1,8
          do i = 1,inlon
            c(i,j,n) = c(i,j,n)*pfrac(i0+i*nint,jj)
          enddo ! i = 1,inlon
        enddo ! n = 1,8
        do i = 1,inlon
          c(i,j,9) = c(i,j,9)*pfrac(i0+i*nint,jj)+ &
            (1._r8-pfrac(i0+i*nint,jj))*c(i,j,9)* &
            (dlatm*real8/(10._r8*dtr))**2
        enddo ! i = 1,inlon
      enddo ! j=1,inlat
    endif ! nint
  end subroutine stenmd
!--------------------------------------------------------------------
  subroutine solver(cofum,c0)
   use edyn_mudmod, only: mudmod
   use edyn_muh2cr, only: muh
!
! Call mudpack to solve PDE. Solution is returned in rim:
!     real,dimension(nmlonp1,nmlat,2) :: rim
! Mudpack solvers:
!   isolve = 0  org. mud v5.      (modified stencils neq direct solution)
!   isolve = 1  muh hybrid solver (no convergence => only as direct solver)
!   isolve = 2  modified mud      (residual calculated with unmodified stencils
!                                  same solution as with direct solver, if same
!                                  coefficient matrix is used)
! Only isolve==2 is supported in edynamo.
!
! Args:
    real(r8),dimension(nmlon0,nmlat0,9),intent(in) :: cofum
    real(r8),intent(in) :: c0(nmlon0,nmlat0,10)
!
! Local:
    integer :: i,j,jntl,ier,isolve
    real(r8) :: l2norm
    real(r8),dimension(nmlon0,nmlat0)   :: ressolv
    real(r8),dimension(nmlon0,nmlat0,9) :: cofum_solv

! Module data:
! real,dimension(nmlonp1,nmlat,2) :: rim_glb ! pde solver output

    jntl = 0
    ier = 0
    isolve = 2
    call mudmod(rim_glb,phisolv,jntl,isolve,res_nlev,ier)
    if (ier < 0 ) then       ! not converged
      if (masterproc) write(iulog,*) 'solver: use muh direct solver'
      call muh(rim_glb,nmlon,nmlat,res_nlev,jntl)
    endif

    l2norm=0._r8
    ressolv = 0.0_r8
    do j = 1,nmlat0
      do i = 1,nmlon0-1
        cofum_solv(i,j,:)=  cofum(i,j,:)
!
! fields: phisolv(0:nmlonp1,0:nmlat+1)       ! 2d solution/ electric potential
!
        ressolv(i,j) = (                    &
        cofum_solv(i,j,1)*phisolv(i+1,j)+   &
        cofum_solv(i,j,2)*phisolv(i+1,j+1)+ &
        cofum_solv(i,j,3)*phisolv(i,j+1)+   &
        cofum_solv(i,j,4)*phisolv(i-1,j+1)+ &
        cofum_solv(i,j,5)*phisolv(i-1,j)+   &
        cofum_solv(i,j,6)*phisolv(i-1,j-1)+ &
        cofum_solv(i,j,7)*phisolv(i,j-1)+   &
        cofum_solv(i,j,8)*phisolv(i+1,j-1)+ &
        cofum_solv(i,j,9)*phisolv(i,j))

        ressolv(i,j) = c0(i,j,10)-ressolv(i,j)
        l2norm = l2norm + ressolv(i,j)*ressolv(i,j)
      enddo
    enddo
!   write(iulog,*) 'L2norm (global root task) = ',l2norm

  end subroutine solver
!--------------------------------------------------------------------
end module edyn_solve
