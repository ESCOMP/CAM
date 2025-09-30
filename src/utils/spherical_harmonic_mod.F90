module spherical_harmonic_mod
!======================================================================
!
! Purpose: Compute and make use of Spherical Harmonic Analysis on physgrid
!
!    This module implements 1 data structures for the spectral analysis
!    and synthesis of spherical harmonic function.
!
!    SphericalHarmonic_t:    For the analysis/synthesis of spherical harmonic
!                            functions and amplitudes on a 2D grid of points 
!                            distributed over the surface of a sphere.
!
!         The SphericalHarmonic_t computes global integrals to compute basis
!         amplitudes. 
!
! USAGE:
!
!    Compute Spherical Harmonic amplitudes and synthesize values on 2D/3D physgrid
!
!         Usage: type(SphericalHarmonic_t):: SH
!         =========================================
!           call SH%init(nmax,nbas)
!           ------------------
!               - Initialize the data structure for the given spherical 
!                 truncation 'nmax' and return 'nbas',the number of spherical 
!                 harmonic basis functions.
!
!               Arguments:
!                 integer ,intent(in ):: nmax    -Number of meridional modes
!                 integer ,intent(out):: nbas    -Total number spherical harmonic functions
!
!           call SH%calc_amps(Gdata,Bamp)
!           -----------------------------
!               - For the initialized SphericalHarmonic_t; Given Gdata() values on 
!                 the physgrid, compute the harmonic basis amplitudes Bamp().
!
!               Interface: 2D data on the physgrid
!                 real(r8),intent(in ):: Gdata(pcols,begchunk:endchunk)
!                 real(r8),intent(out):: Bamp (nbas)
!
!               Interface: 3D data on the physgrid
!                 real(r8),intent(in ):: Gdata(pcols,pver,begchunk:endchunk)
!                 real(r8),intent(out):: Bamp (nbas,pver)
!
!           call SH%eval_grid(Bamp,Gdata)
!           -----------------------------
!               - For the initialized SphericalHarmonic_t; Given Bamp() spherical 
!                 harmonic basis amplitudes, compute the Gdata() values on the physgrid.
!
!               Interface: 2D data on the physgrid
!                 real(r8),intent(in ):: Bamp (nbas)
!                 real(r8),intent(out):: Gdata(pcols,begchunk:endchunk)
!
!               Interface: 3D data on the physgrid
!                 real(r8),intent(in ):: Bamp (nbas,pver)
!                 real(r8),intent(out):: Gdata(pcols,pver,begchunk:endchunk)
!
!======================================================================

  use shr_kind_mod,    only: r8=>SHR_KIND_R8
  use phys_grid,       only: get_ncols_p, get_rlat_p, get_rlon_p, get_wght_all_p, get_nlcols_p
  use ppgrid,          only: begchunk, endchunk, pcols
  use shr_reprosum_mod,only: shr_reprosum_calc
  use cam_abortutils,  only: endrun, handle_allocate_error
  use spmd_utils,      only: mpicom
  use physconst,       only: pi
  use phys_grid,       only: ngcols_p => num_global_phys_cols
  use cam_logfile,     only: iulog

  implicit none
  private

  public :: SphericalHarmonic_t

  ! Type definitions
  !-------------------
  type SphericalHarmonic_t
     private
     integer             :: nmax
     integer             :: nbas
     real(r8),allocatable:: area (:,:)
     real(r8),allocatable:: basis(:,:,:)
    contains
     procedure,pass:: init      => init_SphericalHarmonic
     generic,public:: calc_amps => calc_SphericalHarmonic_2Damps, &
                                   calc_SphericalHarmonic_3Damps
     generic,public:: eval_grid => eval_SphericalHarmonic_2Dgrid, &
                                   eval_SphericalHarmonic_3Dgrid
     procedure,private,pass:: calc_SphericalHarmonic_2Damps
     procedure,private,pass:: calc_SphericalHarmonic_3Damps
     procedure,private,pass:: eval_SphericalHarmonic_2Dgrid
     procedure,private,pass:: eval_SphericalHarmonic_3Dgrid
     procedure, pass :: final => final_SphericalHarmonic
  end type SphericalHarmonic_t

  real(r8), parameter :: halfPI = 0.5_r8*pi
  real(r8), parameter :: twoPI  = 2.0_r8*pi
  real(r8), parameter :: fourPI = 4.0_r8*pi
  real(r8), parameter :: qrtrPI = 0.25_r8*pi
  real(r8), parameter :: invSqrt4pi = 1._r8/sqrt(fourPI)

contains
    !=======================================================================
    subroutine init_SphericalHarmonic(this,I_nmax,O_nbas)
      !
      ! init_SphericalHarmonic: Initialize the SphericalHarmonic data structure 
      !                 for the physics grid. It is assumed that the domain
      !                 of these gridpoints spans the surface of the sphere.
      !                 The representation of basis functions is
      !                 normalized w.r.t integration over the sphere.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_t) :: this
      integer ,intent(in ):: I_nmax
      integer ,intent(out):: O_nbas
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Clons(:,:)
      real(r8),allocatable:: Clats(:,:)
      real(r8),allocatable:: Bcoef(:)
      real(r8),allocatable:: Bsum (:,:)
      real(r8),allocatable:: Bnorm (:)
      real(r8),allocatable:: Bamp  (:)
      real(r8):: area(pcols),rlat,rlon
      real(r8):: Pnm

      integer :: nbas_e,nbas_o
      integer :: mm,nn,n2,nb,lchnk,ncols,cc

      integer :: nlcols, count, astat
      character(len=*), parameter :: subname = 'init_SphericalHarmonic'

      if (I_nmax<1) then
         call endrun('SphericalHarmonic%init: ERROR I_nmax must be greater than 0')
      end if

      ! Allocate space
      !-----------------
      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)

      O_nbas = (I_nmax+1)**2
      this%nmax = I_nmax
      this%nbas = O_nbas
      allocate(this%area (pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'this%area')
      allocate(this%basis(pcols,begchunk:endchunk,O_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%basis')
      this%area (:,:)   = 0._r8
      this%basis(:,:,:) = 0._r8

      nlcols = get_nlcols_p()

      allocate(Clons(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'Clons')
      allocate(Clats(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'Clats')
      allocate(Bcoef(O_nbas/2+1), stat=astat)
      call handle_allocate_error(astat, subname, 'Bcoef')
      allocate(Bsum (nlcols,1), stat=astat)
      call handle_allocate_error(astat, subname, 'Bsum')
      allocate(Bamp (1), stat=astat)
      call handle_allocate_error(astat, subname, 'Bamp')
      allocate(Bnorm(1), stat=astat)
      call handle_allocate_error(astat, subname, 'Bnorm')

      ! Save a copy of the area weights for each ncol gridpoint
      ! and convert Latitudes to SP->NP colatitudes in radians
      !-------------------------------------------------------
      do lchnk=begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        call get_wght_all_p(lchnk, ncols, area)
        do cc = 1,ncols
          rlat=get_rlat_p(lchnk,cc)
          rlon=get_rlon_p(lchnk,cc)
          this%area(cc,lchnk) = area(cc)
          Clons    (cc,lchnk) = rlon
          Clats    (cc,lchnk) = rlat + halfPI
        end do
      end do

      ! Add first basis for the mean values.
      !------------------------------------------
      this%nbas = 1
      this%basis(:,begchunk:endchunk,1) = invSqrt4pi

      ! Loop over the remaining meridional modes
      !------------------------------------------
      do nb=2,(this%nmax+1)
        nn = nb-1

        ! Add the m=0 mode first
        !------------------------
        this%nbas = this%nbas + 1
        call sh_gen_basis_coefs(nn,0,Bcoef)
        do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do cc = 1,ncols
            call sh_create_basis(nn,0,Clats(cc,lchnk),Bcoef,this%basis(cc,lchnk,this%nbas))
          end do
        end do

        ! Now loop over zonal modes mm=1,nn
        ! and add even/odd basis functions
        !----------------------------------
        do mm=1,nn
          nbas_e = this%nbas + 1
          nbas_o = this%nbas + 2
          this%nbas = this%nbas + 2
          call sh_gen_basis_coefs(nn,mm,Bcoef)
          do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols
              call sh_create_basis(nn,mm,Clats(cc,lchnk),Bcoef,Pnm)
              this%basis(cc,lchnk,nbas_e) = Pnm*cos(mm*Clons(cc,lchnk))
              this%basis(cc,lchnk,nbas_o) = Pnm*sin(mm*Clons(cc,lchnk))
            end do
          end do
        end do
      end do ! nn=2,this%nbas
      O_nbas = this%nbas

      !-------------------------------------------------------------
      ! The Discrete basis representation needs to be orthogonalized
      !  Grahm-Schmidt Orthogonization
      !-------------------------------------------------------------

      ! Numerically normalize the gravest functon
      ! as the first orthonormal basis function.
      !-------------------------------------------
      count = 0
      Bsum(:,:) = 0._r8
      do lchnk=begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        do cc = 1,ncols
          count = count+1
          Bsum(count,1) = this%basis(cc,lchnk,1)*this%basis(cc,lchnk,1)*this%area(cc,lchnk)
        end do
      end do

      call shr_reprosum_calc(Bsum, Bnorm, count, nlcols, 1, gbl_count=ngcols_p, commid=mpicom)

      do lchnk=begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        this%basis(:ncols,lchnk,1) = this%basis(:ncols,lchnk,1)/sqrt(Bnorm(1))
      end do

      ! Loop over the remaining basis functions
      !-------------------------------------------
      do nn=2,this%nbas

        ! Remove contributions from exisiting set of orthonormal functions
        !------------------------------------------------------------------
        do n2=1,(nn-1)
          count = 0
          Bsum(:,:) = 0._r8
          do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols
              count = count+1
              Bsum(count,1) = this%basis(cc,lchnk,nn)*this%basis(cc,lchnk,n2)*this%area(cc,lchnk)
            end do
          end do

          call shr_reprosum_calc(Bsum, Bamp, count, nlcols, 1, gbl_count=ngcols_p, commid=mpicom)

          do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            this%basis(:ncols,lchnk,nn) = this%basis(:ncols,lchnk,nn) - Bamp(1)*this%basis(:ncols,lchnk,n2)
          end do
        end do ! n2=1,(nn-1)

        ! Normalize the result for the newest member of the orthonomal set
        !--------------------------------------------------------------------
        count = 0
        Bsum(:,:) = 0._r8
        do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do cc = 1,ncols
            count = count+1
            Bsum(count,1) = this%basis(cc,lchnk,nn)*this%basis(cc,lchnk,nn)*this%area(cc,lchnk)
          end do
        end do
  
        call shr_reprosum_calc(Bsum, Bnorm, count, nlcols, 1, gbl_count=ngcols_p, commid=mpicom)
  
        do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          this%basis(:ncols,lchnk,nn) = this%basis(:ncols,lchnk,nn)/sqrt(Bnorm(1))
        end do
      end do ! nn=2,this%nbas

      !DIAG:  Check for blatent orthogonality errors
      !-----------------------------------------------
!   IF(.TRUE.) THEN
!      write(iulog,*)  'PFC:  ORTHONORM CHECK:'
!      do nn=2,this%nbas
!        do n2=1,(nn-1)
!          count = 0
!          Bsum(:,:) = 0._r8
!          do lchnk=begchunk,endchunk
!            ncols = get_ncols_p(lchnk)
!            do cc = 1,ncols
!              count = count+1
!              Bsum(count,1) = this%basis(cc,lchnk,nn)*this%basis(cc,lchnk,n2)*this%area(cc,lchnk)
!            end do
!          end do
!          call shr_reprosum_calc(Bsum, Bamp, count, nlcols, 1, gbl_count=ngcols_p, commid=mpicom)
!          if(abs(Bamp(1)).gt.1.d-5) then
!            write(iulog,*) 'PFC: *** nn=',nn,' n2=',n2,'  Bamp=',Bamp(1)
!          endif
!        end do
!      end do ! nn=2,this%nbas
!      write(iulog,*)  'PFC:  ORTHONORM CHECK: done'
!   ENDIF

      ! End Routine
      !------------
      deallocate(Clons)
      deallocate(Clats)
      deallocate(Bcoef)
      deallocate(Bsum )
      deallocate(Bamp )
      deallocate(Bnorm)
    end subroutine init_SphericalHarmonic
    !=======================================================================


    !=======================================================================
    subroutine calc_SphericalHarmonic_2Damps(this,I_Gdata,O_Bamp)
      !
      ! calc_SphericalHarmonic_2Damps: Given 2D data values for the ncol gridpoints,
      !                        compute the spherical harmonic basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_t) :: this
      real(r8),intent(in ) :: I_Gdata(pcols,begchunk:endchunk)
      real(r8),intent(out) :: O_Bamp(:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable :: Csum(:,:)
      real(r8),allocatable :: Bamp(:)
      integer :: nn,n2,ncols,lchnk,cc
      integer :: nlcols, count, astat

      character(len=*), parameter :: subname = 'calc_SphericalHarmonic_2Damps'

      nlcols = get_nlcols_p()

      allocate(Bamp(this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bamp')
      allocate(Csum(nlcols, this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Csum')
      Csum(:,:) = 0._r8

      ! Compute Covariance with input data and basis functions
      !--------------------------------------------------------
      do nn= 1,this%nbas
        count = 0
        do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do cc = 1,ncols
            count=count+1
            Csum(count,nn) = I_Gdata(cc,lchnk)*this%basis(cc,lchnk,nn)*this%area(cc,lchnk)
          end do
        end do
      end do

      call shr_reprosum_calc(Csum, Bamp, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

      ! Output the amplitudes
      !--------------------------
      do nn=1,this%nbas
        O_Bamp(nn) = Bamp(nn)
      end do

      ! End Routine
      !------------
      deallocate(Csum)
      deallocate(Bamp)

    end subroutine calc_SphericalHarmonic_2Damps
    !=======================================================================


    !=======================================================================
    subroutine calc_SphericalHarmonic_3Damps(this,I_Gdata,O_Bamp)
      !
      ! calc_SphericalHarmonic_3Damps: Given 3D data values for the ncol,nlev gridpoints,
      !                        compute the zonal mean basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_t) :: this
      real(r8),intent(in ):: I_Gdata(:,:,begchunk:)
      real(r8),intent(out):: O_Bamp (:,:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Csum (:,:)
      real(r8),allocatable:: Bamp (:)
      integer:: nn,n2,ncols,lchnk,cc
      integer:: Nsum,ns,ll
      integer :: nlcols, count, astat

      integer :: nlev
      character(len=*), parameter :: subname = 'calc_SphericalHarmonic_3Damps'

      nlev = size(I_Gdata,dim=2)

      nlcols = get_nlcols_p()
      allocate(Bamp(this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bamp')
      allocate(Csum(nlcols, this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Csum')

      Csum(:,:) = 0._r8
      O_Bamp(:,:) = 0._r8

      ! Compute Covariance with input data and basis functions
      !--------------------------------------------------------
      do ll= 1,nlev

         Csum(:,:) = 0._r8
         Bamp(:) = 0._r8

         do nn= 1,this%nbas
            count = 0
            do lchnk=begchunk,endchunk
               ncols = get_ncols_p(lchnk)
               do cc = 1,ncols
                  count=count+1
                  Csum(count,nn) = I_Gdata(cc,ll,lchnk)*this%basis(cc,lchnk,nn)*this%area(cc,lchnk)
               end do
            end do
         end do

         call shr_reprosum_calc(Csum, Bamp, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

         ! Output the amplitudes
         !--------------------------
         do nn=1,this%nbas
            O_Bamp(nn,ll) = Bamp(nn)
         end do

      end do

      ! End Routine
      !------------
      deallocate(Csum)
      deallocate(Bamp)

    end subroutine calc_SphericalHarmonic_3Damps
    !=======================================================================


    !=======================================================================
    subroutine eval_SphericalHarmonic_2Dgrid(this,I_Bamp,O_Gdata)
      !
      ! eval_SphericalHarmonic_2Dgrid: Given the zonal mean basis amplitudes,
      !                        compute 2D data values for the ncol gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_t) :: this
      real(r8),intent(in ):: I_Bamp (:)
      real(r8),intent(out):: O_Gdata(pcols,begchunk:endchunk)
      !
      ! Local Values
      !--------------
      integer:: nn,ncols,lchnk,cc

      O_Gdata(:,:) = 0._r8

      ! Construct grid values from basis amplitudes.
      !--------------------------------------------------

      do nn=1,this%nbas
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols
               O_Gdata(cc,lchnk) = O_Gdata(cc,lchnk) + (I_Bamp(nn)*this%basis(cc,lchnk,nn))
            end do
         end do
      end do

    end subroutine eval_SphericalHarmonic_2Dgrid
    !=======================================================================


    !=======================================================================
    subroutine eval_SphericalHarmonic_3Dgrid(this,I_Bamp,O_Gdata)
      !
      ! eval_SphericalHarmonic_3Dgrid: Given the zonal mean basis amplitudes,
      !                      compute 3D data values for the ncol,nlev gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_t) :: this
      real(r8),intent(in ):: I_Bamp (:,:)
      real(r8),intent(out):: O_Gdata(:,:,begchunk:)
      !
      ! Local Values
      !--------------
      integer:: nn,ncols,lchnk,cc
      integer:: ll

      integer :: nlev
      nlev = size(O_Gdata,dim=2)

      O_Gdata(:,:,:) = 0._r8

      ! Construct grid values from basis amplitudes.
      !--------------------------------------------------

      do ll = 1,nlev
         do nn=1,this%nbas
            do lchnk=begchunk,endchunk
               ncols = get_ncols_p(lchnk)
               do cc = 1,ncols
                  O_Gdata(cc,ll,lchnk) = O_Gdata(cc,ll,lchnk) + (I_Bamp(nn,ll)*this%basis(cc,lchnk,nn))
               end do
            end do
         end do
      end do

    end subroutine eval_SphericalHarmonic_3Dgrid
    !=======================================================================


    !=======================================================================
    subroutine final_SphericalHarmonic(this)
      class(SphericalHarmonic_t) :: this

      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)

    end subroutine final_SphericalHarmonic
    !=======================================================================


    !=======================================================================
    subroutine Invert_Matrix(I_Mat,Nbas,O_InvMat)
      !
      ! Invert_Matrix: Given the NbasxNbas matrix, calculate and return
      !                the inverse of the matrix.
      !
      ! Implemented with the LAPACK DGESV routine.
      !
      !====================================================================
      !
      ! Passed Variables
      !------------------
      real(r8), intent(inout) :: I_Mat(:,:)  ! input matrix contains P*L*U
                                             ! decomposition on output
      integer,  intent(in)    :: Nbas
      real(r8), intent(out)   :: O_InvMat(:,:)
      !
      ! Local Values
      !-------------
      integer, allocatable :: Indx(:)  ! pivot indices
      integer :: astat, ii
      character(len=*), parameter :: subname = 'Invert_Matrix'
      character(len=80) :: msg

      external DGESV

      ! Allocate work space
      !---------------------
      allocate(Indx(Nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Indx')

      ! Initialize the inverse array with the identity matrix
      !-------------------------------------------------------
      O_InvMat(:,:) = 0._r8
      do ii=1,Nbas
        O_InvMat(ii,ii) = 1._r8
      end do

      call DGESV(Nbas, Nbas, I_Mat, Nbas, Indx, O_InvMat, Nbas, astat)

      if (astat < 0) then
         write(msg, '(a, i1, a)') 'argument # ', abs(astat), ' has an illegal value'
         call endrun(subname//': DGESV error return: '//msg)
      else if (astat > 0) then
         call endrun(subname//': DGESV error return: matrix is singular')
      end if

      deallocate(Indx)

    end subroutine Invert_Matrix
    !=======================================================================

    !=======================================================================
    ! legacy spherepack routines
    !=======================================================================
    subroutine sh_gen_basis_coefs(nn,mm,cp)
      !
      ! spherepack alfk
      !
      ! dimension of           real cp(nn/2 + 1)
      ! arguments
      !
      ! purpose                computes fourier coefficients in the trigonometric series
      !                        representation of the normalized associated
      !                        legendre function pbar(nn,mm,theta) for use by
      !                        sh_gen_basis_coefs in calculating pbar(nn,mm,theta).
      !
      !                        first define the normalized associated
      !                        legendre functions
      !
      !                        pbar(mm,nn,theta) = sqrt((2*nn+1)*factorial(nn-mm)
      !                        /(2*factorial(nn+mm)))*sin(theta)**mm/(2**nn*
      !                        factorial(nn)) times the (nn+mm)th derivative of
      !                        (x**2-1)**nn with respect to x=cos(theta)
      !
      !                        where theta is colatitude.
      !
      !                        then subroutine sh_gen_basis_coefs computes the coefficients
      !                        cp(k) in the following trigonometric
      !                        expansion of pbar(m,n,theta).
      !
      !                        1) for n even and m even, pbar(mm,nn,theta) =
      !                           .5*cp(1) plus the sum from k=1 to k=nn/2
      !                           of cp(k+1)*cos(2*k*th)
      !
      !                        2) for nn even and mm odd, pbar(mm,nn,theta) =
      !                           the sum from k=1 to k=nn/2 of
      !                           cp(k)*sin(2*k*th)
      !
      !                        3) for n odd and m even, pbar(mm,nn,theta) =
      !                           the sum from k=1 to k=(nn+1)/2 of
      !                           cp(k)*cos((2*k-1)*th)
      !
      !                        4) for nn odd and mm odd,  pbar(mm,nn,theta) =
      !                           the sum from k=1 to k=(nn+1)/2 of
      !                           cp(k)*sin((2*k-1)*th)
      !
      ! arguments
      !
      ! on input               nn
      !                          nonnegative integer specifying the degree of
      !                          pbar(nn,mm,theta)
      !
      !                        mm
      !                          is the order of pbar(nn,mm,theta). mm can be
      !                          any integer however cp is computed such that
      !                          pbar(nn,mm,theta) = 0 if abs(m) is greater
      !                          than nn and pbar(nn,mm,theta) = (-1)**mm*
      !                          pbar(nn,-mm,theta) for negative mm.
      !
      ! on output              cp
      !                          array of length (nn/2)+1
      !                          which contains the fourier coefficients in
      !                          the trigonometric series representation of
      !                          pbar(nn,mm,theta)
      !
      ! special conditions     none
      !
      ! algorithm              the highest order coefficient is determined in
      !                        closed form and the remainig coefficients are
      !                        determined as the solution of a backward
      !                        recurrence relation.
      !
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      integer ,intent(in ):: nn
      integer ,intent(in ):: mm
      real(r8),intent(out):: cp(nn/2+1)
      !
      ! Local Values
      !----------------
      real(r8):: fnum,fnmh
      real(r8):: pm1
      real(r8):: t1,t2
      real(r8):: fden
      real(r8):: cp2
      real(r8):: fnnp1
      real(r8):: fnmsq
      real(r8):: fk
      real(r8):: a1,b1,C1
      integer :: ma,nmms2,nex
      integer :: ii,jj

      real(r8),parameter:: SC10=1024._r8
      real(r8),parameter:: SC20=SC10*SC10
      real(r8),parameter:: SC40=SC20*SC20

      cp(1) = 0._r8
      ma = abs(mm)
      if(ma>nn) return

      if((nn-1)<0) then
        cp(1) = sqrt(2._r8)
        return
      elseif((nn-1)==0) then
        if(ma/=0) then
          cp(1) = sqrt(.75_r8)
          if(mm==-1) cp(1) = -cp(1)
        else
          cp(1) = sqrt(1.5_r8)
        endif
        return
      else
        if(mod(nn+ma,2)/=0) then
          nmms2 = (nn-ma-1)/2
          fnum  = nn + ma + 2
          fnmh  = nn - ma + 2
          pm1   = -1._r8
        else
          nmms2 = (nn-ma)/2
          fnum  = nn + ma + 1
          fnmh  = nn - ma + 1
          pm1   = 1._r8
        endif
      endif

      t1   = 1._r8/SC20
      nex  = 20
      fden = 2._r8
      if(nmms2>=1) then
        do ii = 1,nmms2
          t1 = fnum*t1/fden
          if (t1>SC20) then
            t1  = t1/SC40
            nex = nex + 40
          endif
          fnum = fnum + 2._r8
          fden = fden + 2._r8
        end do
      endif

      if(mod(ma/2,2)/=0) then
        t1 = -t1/2._r8**(nn-1-nex)
      else
        t1 =  t1/2._r8**(nn-1-nex)
      endif
      t2 = 1._r8
      if(ma/=0) then
        do ii = 1,ma
          t2   = fnmh*t2/ (fnmh+pm1)
          fnmh = fnmh + 2._r8
        end do
      endif

      cp2   = t1*sqrt((nn+.5_r8)*t2)
      fnnp1 = nn*(nn+1)
      fnmsq = fnnp1 - 2._r8*ma*ma

      if((mod(nn,2)==0).and.(mod(ma,2)==0)) then
        jj = 1+(nn+1)/2
      else
        jj = (nn+1)/2
      endif

      cp(jj) = cp2
      if(mm<0) then
        if(mod(ma,2)/=0) cp(jj) = -cp(jj)
      endif
      if(jj<=1) return

      fk = nn
      a1 = (fk-2._r8)*(fk-1._r8) - fnnp1
      b1 = 2._r8* (fk*fk-fnmsq)
      cp(jj-1) = b1*cp(jj)/a1

      jj = jj - 1
      do while(jj>1)
        fk = fk - 2._r8
        a1 = (fk-2._r8)*(fk-1._r8) - fnnp1
        b1 = -2._r8*(fk*fk-fnmsq)
        c1 = (fk+1._r8)*(fk+2._r8) - fnnp1
        cp(jj-1) = -(b1*cp(jj)+c1*cp(jj+1))/a1
        jj = jj - 1
     end do

    end subroutine sh_gen_basis_coefs
    !=======================================================================

    !=======================================================================
    subroutine sh_create_basis(nn,mm,theta,cp,pb)
      !
      ! spherepack lfpt
      !
      ! dimension of
      ! arguments
      !                        cp((nn/2)+1)
      !
      ! purpose                routine sh_create_basis uses coefficients computed by
      !                        routine sh_gen_basis_coefs to compute the
      !                        normalized associated legendre function pbar(nn,mm,theta)
      !                        at colatitude theta.
      !
      ! arguments
      !
      ! on input               nn
      !                          nonnegative integer specifying the degree of
      !                          pbar(nn,mm,theta)
      !                        mm
      !                          is the order of pbar(nn,mm,theta). mm can be
      !                          any integer however pbar(nn,mm,theta) = 0
      !                          if abs(mm) is greater than nn and
      !                          pbar(nn,mm,theta) = (-1)**mm*pbar(nn,-mm,theta)
      !                          for negative mm.
      !
      !                        theta
      !                          colatitude in radians
      !
      !                        cp
      !                          array of length (nn/2)+1
      !                          containing coefficients computed by routine
      !                          sh_gen_basis_coefs
      !
      ! on output              pb
      !                          variable containing pbar(n,m,theta)
      !
      ! special conditions     calls to routine sh_create_basis must be preceded by an
      !                        appropriate call to routine sh_gen_basis_coefs.
      !
      ! algorithm              the trigonometric series formula used by
      !                        routine sh_create_basis to calculate pbar(nn,mm,theta) at
      !                        colatitude theta depends on mm and nn as follows:
      !
      !                           1) for nn even and mm even, the formula is
      !                              .5*cp(1) plus the sum from k=1 to k=n/2
      !                              of cp(k)*cos(2*k*theta)
      !                           2) for nn even and mm odd. the formula is
      !                              the sum from k=1 to k=nn/2 of
      !                              cp(k)*sin(2*k*theta)
      !                           3) for nn odd and mm even, the formula is
      !                              the sum from k=1 to k=(nn+1)/2 of
      !                              cp(k)*cos((2*k-1)*theta)
      !                           4) for nn odd and mm odd, the formula is
      !                              the sum from k=1 to k=(nn+1)/2 of
      !                              cp(k)*sin((2*k-1)*theta)
      !
      !=====================================================================
      integer, intent(in) :: nn,mm
      real(r8), intent(in) :: theta
      real(r8), intent(in) :: cp(:)
      real(r8), intent(out) :: pb

      real(r8) :: cdt
      real(r8) :: sdt
      real(r8) :: ct
      real(r8) :: st
      real(r8) :: summ
      real(r8) :: cth

      integer:: ma,nmod,mmod,kdo
      integer:: kp1,kk

      pb = 0._r8
      ma = abs(mm)
      if(ma>nn) return

      if(nn<=0) then
        if(ma<=0) then
          pb = sqrt(.5_r8)
          return
        endif
      endif

      nmod = mod(nn,2)
      mmod = mod(ma,2)

      if(nmod<=0) then
        if(mmod<=0) then
          kdo = nn/2 + 1
          cdt = cos(theta+theta)
          sdt = sin(theta+theta)
          ct  = 1._r8
          st  = 0._r8
          summ = .5_r8*cp(1)
          do kp1 = 2,kdo
            cth = cdt*ct - sdt*st
            st  = sdt*ct + cdt*st
            ct  = cth
            summ = summ + cp(kp1)*ct
          end do
          pb = summ
          return
        endif
        kdo = nn/2
        cdt = cos(theta+theta)
        sdt = sin(theta+theta)
        ct  = 1._r8
        st  = 0._r8
        summ = 0._r8
        do kk = 1,kdo
          cth = cdt*ct - sdt*st
          st  = sdt*ct + cdt*st
          ct  = cth
          summ = summ + cp(kk)*st
        end do
        pb = summ
        return
      endif

      kdo = (nn+1)/2
      if(mmod<=0) then
        cdt =  cos(theta+theta)
        sdt =  sin(theta+theta)
        ct  =  cos(theta)
        st  = -sin(theta)
        summ = 0._r8
        do kk = 1,kdo
          cth = cdt*ct - sdt*st
          st  = sdt*ct + cdt*st
          ct  = cth
          summ = summ + cp(kk)*ct
        end do
        pb = summ
        return
      endif

      cdt =  cos(theta+theta)
      sdt =  sin(theta+theta)
      ct  =  cos(theta)
      st  = -sin(theta)
      summ = 0._r8
      do kk = 1,kdo
        cth = cdt*ct - sdt*st
        st  = sdt*ct + cdt*st
        ct  = cth
        summ = summ + cp(kk)*st
      end do
      pb = summ

    end subroutine sh_create_basis
    !=======================================================================

    !=======================================================================
    subroutine sh_create_gaus_grid(nlat,theta,wts,ierr)
      !
      ! spherepack gaqd
      !  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      !  .                                                             .
      !  .                  copyright (c) 2001 by ucar                 .
      !  .                                                             .
      !  .       university corporation for atmospheric research       .
      !  .                                                             .
      !  .                      all rights reserved                    .
      !  .                                                             .
      !  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      !
      !                             February 2002
      !
      !     gauss points and weights are computed using the fourier-newton
      !     described in "on computing the points and weights for
      !     gauss-legendre quadrature", paul n. swarztrauber, siam journal
      !     on scientific computing (DOI 10.1137/S1064827500379690).
      !     This routine is faster and more accurate than older program
      !     with the same name.
      !
      !     computes the nlat gaussian colatitudes and weights
      !     in double precision. the colatitudes are in radians and lie in the
      !     in the interval (0,pi).
      !
      !     input parameters
      !
      !     nlat    the number of gaussian colatitudes in the interval (0,pi)
      !             (between the two poles).  nlat must be greater than zero.
      !
      !     output parameters
      !
      !     theta   a double precision array with length nlat
      !             containing the gaussian colatitudes in
      !             increasing radians on the interval (0,pi).
      !
      !     wts     a double precision array with lenght nlat
      !             containing the gaussian weights.
      !
      !     ierror = 0 no errors
      !            = 1 if nlat<=0
      !
      !===================================================================
      !
      ! Passed variables
      !-----------------
      integer ,intent(in ) :: nlat
      real(r8),intent(out) :: theta(nlat)
      real(r8),intent(out) :: wts(nlat)
      integer ,intent(out) :: ierr
      !
      ! Local Values
      !-------------
      real(r8):: sgnd
      real(r8):: xx,dtheta,dthalf
      real(r8):: cmax,zprev,zlast,zero,zhold,pb,dpb,dcor,summ,cz
      integer :: mnlat,ns2,nhalf,nix,it,ii

      real(r8), parameter :: eps = epsilon(1._r8)

      ! check work space length
      !------------------------
      if(nlat<=0) then
        ierr = 1
        return
      endif
      ierr = 0

      ! compute weights and points analytically when nlat=1,2
      !-------------------------------------------------------
      if(nlat==1) then
        theta(1) = acos(0._r8)
        wts  (1) = 2._r8
        return
      elseif(nlat==2) then
        xx       = sqrt(1._r8/3._r8)
        theta(1) = acos( xx)
        theta(2) = acos(-xx)
        wts  (1) = 1._r8
        wts  (2) = 1._r8
        return
      endif

      ! Proceed for nlat > 2
      !----------------------
      mnlat = mod(nlat,2)
      ns2   = nlat/2
      nhalf = (nlat+1)/2

      call sh_fourier_coefs_dp(nlat,cz,theta(ns2+1),wts(ns2+1))

      dtheta = halfPI/nhalf
      dthalf = dtheta/2._r8
      cmax   = .2_r8*dtheta

      ! estimate first point next to theta = pi/2
      !-------------------------------------------
      if(mnlat/=0) then
        zero  = halfPI - dtheta
        zprev = halfPI
        nix   = nhalf - 1
      else
        zero = halfPI - dthalf
        nix  = nhalf
      endif

      do while(nix/=0)
         dcor = huge(1._r8)
         it = 0
         do while (abs(dcor) > eps*abs(zero))
            it = it + 1
            ! newton iterations
            !-----------------------
            call sh_legp_dlegp_theta(nlat,zero,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
            dcor = pb/dpb
            if(dcor.ne.0._r8) then
               sgnd = dcor/abs(dcor)
            else
               sgnd = 1._r8
            endif
            dcor = sgnd*min(abs(dcor),cmax)
            zero = zero - dcor
         end do

         theta(nix) = zero
         zhold      = zero

         !   wts(nix) = (nlat+nlat+1)/(dpb*dpb)
         ! yakimiw's formula permits using old pb and dpb
         !--------------------------------------------------
         wts(nix) = (nlat+nlat+1)/ (dpb+pb*dcos(zlast)/dsin(zlast))**2
         nix      = nix - 1
         if(nix==nhalf-1) zero = 3._r8*zero - pi
         if(nix<nhalf-1) zero = zero + zero - zprev
         zprev = zhold
      end do

      ! extend points and weights via symmetries
      !-------------------------------------------
      if(mnlat/=0) then
        theta(nhalf) = halfPI
        call sh_legp_dlegp_theta(nlat,halfPI,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
        wts(nhalf) = (nlat+nlat+1)/ (dpb*dpb)
      endif

      do ii = 1,ns2
        wts  (nlat-ii+1) = wts(ii)
        theta(nlat-ii+1) = pi - theta(ii)
      end do

      summ = 0._r8
      do ii = 1,nlat
        summ = summ + wts(ii)
      end do
      do ii = 1,nlat
        wts(ii) = 2._r8*wts(ii)/summ
      end do

    end subroutine sh_create_gaus_grid
    !=======================================================================

    !=======================================================================
    subroutine sh_fourier_coefs_dp(nn,cz,cp,dcp)
      !
      ! spherepack cpdp
      !
      !     computes the fourier coefficients of the legendre
      !     polynomial p_n^0 and its derivative.
      !     nn is the degree and nn/2 or (nn+1)/2
      !     coefficients are returned in cp depending on whether
      !     nn is even or odd. The same number of coefficients
      !     are returned in dcp. For nn even the constant
      !     coefficient is returned in cz.
      !=====================================================================
      !
      ! Passed variables
      !-----------------
      integer, intent(in) :: nn
      real(r8), intent(out) :: cz
      real(r8), intent(out) :: cp(nn/2+1)
      real(r8), intent(out) :: dcp(nn/2+1)
      !
      ! Local Values
      !--------------
      real(r8):: t1,t2,t3,t4
      integer :: ncp,jj

      ncp = (nn+1)/2
      t1  = -1._r8
      t2  = nn + 1._r8
      t3  = 0._r8
      t4  = nn + nn + 1._r8
      if(mod(nn,2)==0) then
        cp(ncp) = 1._r8
        do jj = ncp,2,-1
          t1 = t1 + 2._r8
          t2 = t2 - 1._r8
          t3 = t3 + 1._r8
          t4 = t4 - 2._r8
          cp(jj-1) = (t1*t2)/ (t3*t4)*cp(jj)
        end do
        t1 = t1 + 2._r8
        t2 = t2 - 1._r8
        t3 = t3 + 1._r8
        t4 = t4 - 2._r8
        cz = (t1*t2)/ (t3*t4)*cp(1)
        do jj = 1,ncp
          dcp(jj) = (jj+jj)*cp(jj)
        end do
      else
        cp(ncp) = 1._r8
        do jj = ncp - 1,1,-1
          t1 = t1 + 2._r8
          t2 = t2 - 1._r8
          t3 = t3 + 1._r8
          t4 = t4 - 2._r8
          cp(jj) = (t1*t2)/ (t3*t4)*cp(jj+1)
        end do
        do jj = 1,ncp
          dcp(jj) = (jj+jj-1)*cp(jj)
        end do
      endif

    end subroutine sh_fourier_coefs_dp
    !=======================================================================

    !=======================================================================
    subroutine sh_legp_dlegp_theta(nn,theta,cz,cp,dcp,pb,dpb)
      !
      ! spherepack tpdp
      !
      !     computes pn(theta) and its derivative dpb(theta) with
      !     respect to theta
      !=====================================================================
      !
      ! Passed variables
      !------------------
      integer, intent(in) :: nn
      real(r8), intent(in) :: theta
      real(r8), intent(in) :: cz
      real(r8), intent(in) :: cp (nn/2+1)
      real(r8), intent(in) :: dcp(nn/2+1)
      real(r8), intent(out) :: pb
      real(r8), intent(out) :: dpb
      !
      ! Local Values
      !--------------
      real(r8):: cdt,sdt,cth,sth,chh
      integer :: kdo,kk

      cdt = dcos(theta+theta)
      sdt = dsin(theta+theta)
      if(mod(nn,2)==0) then
        ! n even
        !----------
        kdo = nn/2
        pb  = .5_r8*cz
        dpb = 0._r8
        if(nn>0) then
          cth = cdt
          sth = sdt
          do kk = 1,kdo
            pb  = pb  +  cp(kk)*cth
            dpb = dpb - dcp(kk)*sth
            chh = cdt*cth - sdt*sth
            sth = sdt*cth + cdt*sth
            cth = chh
          end do
        endif
      else
        ! n odd
        !-----------
        kdo = (nn+1)/2
        pb  = 0._r8
        dpb = 0._r8
        cth = dcos(theta)
        sth = dsin(theta)
        do kk = 1,kdo
          pb  = pb  +  cp(kk)*cth
          dpb = dpb - dcp(kk)*sth
          chh = cdt*cth - sdt*sth
          sth = sdt*cth + cdt*sth
          cth = chh
        end do
      endif

    end subroutine sh_legp_dlegp_theta
    !=======================================================================

end module spherical_harmonic_mod
