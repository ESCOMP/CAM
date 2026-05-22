module ug_spectralmethods_mod
!======================================================================
!
! Purpose: Implement Spherical Harmonic Analysis/Synthesis methods 
!          for Unstructured Grids(UG).
!
!    Spherical Harmonic Analysis/Synthesis:
!    ======================================
!    This module implements 2 classes for the spectral analysis and 
!    synthesis of spherical harmonic functions. For each, the global
!    integrals to compute basis amlitudes are implemented with an area 
!    weighted sum of the UG gridpoints on the surface of the sphere.
!
!    One computes ampitudes via a least-squares fit in the typical mannor. 
!    This approach distributes numerical errors in the basis representation 
!    based on covariances, but it can be very costly for the MPI environment. 
!    The other approach orthogonalizes the spherical harmonic basis functions 
!    via the Gram-Schmidt method. While this is much less computationally 
!    expensive, it propagates numerical errors in the basis representations
!    from the gravest scales down to the higher resolution modes.
!
!      SphericalHarmonic_COV_t:  For the analysis/synthesis of spherical harmonic
!                                functions and amplitudes on a 2D/3D grids of points 
!                                distributed over the surface of a sphere.
!                                Amplitudes are computed based on COV values between 
!                                basis functions.
!
!      SphericalHarmonic_GS_t:   For the analysis/synthesis of spherical harmonic
!                                functions and amplitudes on a 2D/3D grids of points 
!                                distributed over the surface of a sphere.
!                                Orthonormal basis functions are constructed via
!                                Gram-Schmidt. 
!
!    Zonal Mean Analysis/Synthesis:
!    ======================================
!    This module implements 3 classes for the spectral analysis and 
!    synthesis of zonal mean values based on m=0 spherical harmonics.
!
!      ZonalMean_t:    For the analysis/synthesis of zonal mean values
!                      on a 2D/3D grid of points distributed over the
!                      surface of a sphere.
!      ZonalProfile_t: For the analysis/synthesis of zonal mean values
!                      on a meridional grid that spans the latitudes
!                      from SP to NP
!      ZonalAverage_t: To calculate zonal mean values via a simple
!                      area weighted bin-averaging of 2D/3D grid points
!                      assigned to each latitude band.
!
!    NOTE: The weighting of the Zonal Profiles values is scaled such
!          that ZonalMean_t amplitudes can be used to evaluate values
!          on the ZonalProfile_t grid and vice-versa.
!
!          The ZonalMean_t computes global integrals to compute basis
!          amplitudes. For distributed environments the cost of these
!          can be reduced using the The ZonalAverage_t data structures.
!
!    ======================
!    BASIS SAMPLING ERRORS:
!    ======================
!      The variable values on the unstructured gridpoints represent the 
!      average value over the associated gridpoint's area, but the basis
!      values used at these gridpoints are just the spherical harmonics
!      evaluated at that point. These sampling errors can be significant
!      particularly in the polar regions where the range of lat/lon values 
!      can vary greatly over the domain of a finite area. 
!
!      To suppress these errors, during the initailization process, each 
!      spherical harmonic basis functions is averaged over a local circular 
!      domain centered at each gridpoint. The radius of this domain is set 
!      to match the area associated with the given gridpoint. The number of
!      equal-area gridpoints within this circular domain is specified by the 
!      user via an optional argument to the initialization routine. While the
!      representation of the area is approximate, the goal here is for the 
!      gridpoint values of each basis function to also represent the local 
!      area average over the associated gridpoint's area.
!   
!      The gridpoints in the circular domain are constructed in the the form 
!      of concentric rings such that each point has equal area weighting. The
!      users controls the density of samping points by setting the number of 
!      rings. 
!
!         There is always 1 point at center, each additional 
!         ring (kk) contains 8*(kk-1) sample points. The total 
!         number of sample points (Nsamp) is thus:
!
!             SAMPLE_NRING=1   Nsamp=1      (DEFAULT SampleGrid NOT used)
!             SAMPLE_NRING=2   Nsamp=9  
!             SAMPLE_NRING=3   Nsamp=25  
!             SAMPLE_NRING=4   Nsamp=49  
!
!      The SampleGrid is only used during initialization, after that there 
!      is no additional computational cost.
!    
!
! USAGE:
!
!    Spherical Harmonic Analysis/Synthesis:
!    ======================================
!    Compute Spherical Harmonic amplitudes and synthesize values on 2D/3D physgrid.
!    All of the differences between these classes are internal, so methods/interfaces 
!    are the same for the user:
!
!         Usage: type(SphericalHarmonic_COV_t):: SH
!         Usage: type(SphericalHarmonic_GS_t ):: SH
!         =========================================
!           call SH%init(nmax,nbas,SAMPLE_NRING)
!           ------------------
!               - Initialize the data structure for the given spherical 
!                 truncation 'nmax' and return 'nbas',the number of spherical 
!                 harmonic basis functions.
!
!               Arguments:
!                 integer,intent(in )         :: nmax         -Number of meridional modes
!                 integer,intent(out)         :: nbas         -Total number spherical harmonic functions
!                 integer,intent(in ),optional:: SAMPLE_NRING  (Default=1)
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
! USAGE:
!
!    Zonal Mean Analysis/Synthesis:
!    ======================================
!    (1) Compute Zonal mean amplitudes and synthesize values on 2D/3D physgrid
!
!         Usage: type(ZonalMean_t):: ZM
!         =========================================
!           call ZM%init(nbas,SAMPLE_NRING)
!           ---------------------------------
!               - Initialize the data structure with 'nbas' basis functions
!                 for the given physgrid latitudes and areas.
!
!               Arguments:
!                 integer,intent(in)         :: nbas         -Number of m=0 spherical harmonics
!                 integer,intent(in),optional:: SAMPLE_NRING  (Default=1)
!
!           call ZM%calc_amps(Gdata,Bamp)
!           -----------------------------
!               - For the initialized ZonalMean_t; Given Gdata() values on the physgrid,
!                 compute the zonal mean basis amplitudes Bamp().
!
!               Interface: 2D data on the physgrid
!                 real(r8),intent(in ):: Gdata(pcols,begchunk:endchunk)
!                 real(r8),intent(out):: Bamp (nbas)
!
!               Interface: 3D data on the physgrid
!                 real(r8),intent(in ):: Gdata(pcols,pver,begchunk:endchunk)
!                 real(r8),intent(out):: Bamp (nbas,pver)
!
!           call ZM%eval_grid(Bamp,Gdata)
!           -----------------------------
!               - For the initialized ZonalMean_t; Given Bamp() zonal mean basis
!                 amplitudes, compute the Gdata() values on the physgrid.
!
!               Interface: 2D data on the physgrid
!                 real(r8),intent(in ):: Bamp (nbas)
!                 real(r8),intent(out):: Gdata(pcols,begchunk:endchunk)
!
!               Interface: 3D data on the physgrid
!                 real(r8),intent(in ):: Bamp (nbas,pver)
!                 real(r8),intent(out):: Gdata(pcols,pver,begchunk:endchunk)
!
!
!    (2) Compute Zonal mean amplitudes and synthesize values on Zonal profile grid
!
!         Usage: type(ZonalProfile_t):: ZP
!         =========================================
!           call ZP%init(lats,area,nlat,nbas,GEN_GAUSSLATS=.true.)
!           ------------------------------------------------------
!               - Initialize the data structure for the given number of
!                 latitudes. Either use the given Latitudes and weights,
!                 or OPTIONALLY create profile gridpoints and associated
!                 area weights from SP to NP. Then initialize 'nbas' basis
!                 functions for the profile gridpoints.
!                 If the user supplies the lats/area values, the area values must
!                 be correctly scaled such that the global area adds up to 4PI.
!                 Otherwise, the ampitudes between ZonalProfile_t and ZonalMean_t
!                 are not interchangable.
!
!               Arguments:
!                 real(r8),intent(inout):: lats(:) - Latitudes of meridional grid.
!                 real(r8),intent(inout):: area(:) - Area of each meridional gridpoint.
!                 integer ,intent(in)   :: nlat    - Number of meridional gridpoints.
!                 integer ,intent(in)   :: nbas    - Number of m=0 spherical harmonics
!                 logical ,intent(in),optional:: GEN_GAUSLATS - Flag to generate
!                                                               lats/areas values.
!
!           call ZP%calc_amps(Zdata,Bamp)
!           -----------------------------
!               - Given Zdata() on the Zonal profile grid, compute the
!                 zonal basis amplitudes Bamp().
!
!               Interface: 1D data on (nlat) grid
!                 real(r8),intent(in ):: Zdata(nlat) - Meridional Profile data
!                 real(r8),intent(out):: Bamp (nbas) - Zonal Basis Amplitudes
!
!               Interface: 2D data on (nlat,pver) grid
!                 real(r8),intent(in ):: Zdata(nlat,pver) - Meridional Profile data
!                 real(r8),intent(out):: Bamp (nbas,pver) - Zonal Basis Amplitudes
!
!           call ZP%eval_grid(Bamp,Zdata)
!           -----------------------------
!               - Given Bamp() zonal basis amplitudes, evaluate the Zdata()
!                 values on the Zonal profile grid.
!
!               Interface: 1D data on (nlat) grid
!                 real(r8),intent(in ):: Bamp (nbas) - Zonal Basis Amplitudes
!                 real(r8),intent(out):: Zdata(nlat) - Meridional Profile data
!
!               Interface: 2D data on (nlat,pver) grid
!                 real(r8),intent(in ):: Bamp (nbas,pver) - Zonal Basis Amplitudes
!                 real(r8),intent(out):: Zdata(nlat,pver) - Meridional Profile data
!
!    (3) Compute Zonal mean averages (FASTER/LESS-ACCURATE) on Zonal profile grid
!        (For the created zonal profile, just bin average area weighted
!         2D/3D physgrid grid values)
!
!         Usage: type(ZonalAverage_t):: ZA
!         =========================================
!           call ZA%init(lats,area,nlat,GEN_GAUSSLATS=.true.,USE_LINEARWGTS=.true.)
!           --------------------------------------------------
!               - Given the latitude/area for the nlat meridional gridpoints, initialize
!                 the ZonalAverage data structure for computing bin-averaging of physgrid
!                 values. It is assumed that the domain of these gridpoints of the
!                 profile span latitudes from SP to NP.
!
!                 The optional GEN_GAUSSLATS flag allows for the generation of Gaussian
!                 latitude gridpoints. The generated grid over-writes the given values
!                 lats and area passed by the user.
!
!                 The optional USE_LINEARWGTS flag allows for the linear weighting of 
!                 gridpoint values as they are distributed to latitude bins, or to assign 
!                 grid point values to a single bin with equal weighting. 
!                 (default=TRUE for linear weighting)
!
!               Arguments:
!                 real(r8),intent(inout)     :: lats(nlat)      - Latitudes of meridional grid.
!                 real(r8),intent(inout)     :: area(nlat)      - Area of meridional gridpoints.
!                 integer,intent(in)         :: nlat            - Number of meridional gridpoints
!                 logical,intent(in),optional:: GEN_GAUSLATS    - Flag to generate lats/areas values.
!                 logical,intent(in),optional:: USE_LINEARWGTS  - Flag for weighting option.
!
!           call ZA%binAvg(Gdata,Zdata)
!           ---------------------------
!               - For the initialized ZonalAverage_t; Given Gdata() on the physgrid,
!                 compute bin averages and return Zdata() on the Zonal profile grid.
!
!               Interface: 2D data on the physgrid
!                 real(r8),intent(in ):: Gdata(pcols,begchunk:endchunk)
!                 real(r8),intent(out):: Zdata(nlat)
!
!               Interface: 3D data on the physgrid
!                 real(r8),intent(in ):: Gdata(pcols,pver,begchunk:endchunk)
!                 real(r8),intent(out):: Zdata(nlat,pver)
!
!           call ZA%setVal(Zdata,Gdata)
!           ---------------------------
!               - For the initialized ZonalAverage_t; Given Zdata for the profile bins,
!                 set the 2D physgrid values in each bin to the corresponding value.
!
!               Interface: 2D data on the physgrid
!                 real(r8),intent(in ):: Zdata(nlat)
!                 real(r8),intent(out):: Gdata(pcols,begchunk:endchunk)
!
!               Interface: 3D data on the physgrid
!                 real(r8),intent(in ):: Zdata(nlat,pver)
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

  public :: SphericalHarmonic_COV_t
  public :: SphericalHarmonic_GS_t
  public :: ZonalMean_t
  public :: ZonalProfile_t
  public :: ZonalAverage_t

  ! Type definitions
  !-------------------
  type SphericalHarmonic_COV_t
     private
     integer             :: nmax
     integer             :: nbas
     integer             :: nring
     real(r8),allocatable:: area (:,:)
     real(r8),allocatable:: basis(:,:,:)
     real(r8),allocatable:: map  (:,:)
     integer ,allocatable:: bindx(:,:)
    contains
     procedure,pass:: init      => init_SphericalHarmonic_COV
     procedure,pass:: final     => final_SphericalHarmonic_COV
     generic,public:: calc_amps => calc_SphericalHarmonic_2Damps_COV, &
                                   calc_SphericalHarmonic_3Damps_COV
     generic,public:: eval_grid => eval_SphericalHarmonic_2Dgrid_COV, &
                                   eval_SphericalHarmonic_3Dgrid_COV
     procedure,private,pass:: calc_SphericalHarmonic_2Damps_COV
     procedure,private,pass:: calc_SphericalHarmonic_3Damps_COV
     procedure,private,pass:: eval_SphericalHarmonic_2Dgrid_COV
     procedure,private,pass:: eval_SphericalHarmonic_3Dgrid_COV
  end type SphericalHarmonic_COV_t

  type SphericalHarmonic_GS_t
     private
     integer             :: nmax
     integer             :: nbas
     integer             :: nring
     real(r8),allocatable:: area (:,:)
     real(r8),allocatable:: basis(:,:,:)
    contains
     procedure,pass:: init      => init_SphericalHarmonic_GS
     procedure,pass:: final     => final_SphericalHarmonic_GS
     generic,public:: calc_amps => calc_SphericalHarmonic_2Damps_GS, &
                                   calc_SphericalHarmonic_3Damps_GS
     generic,public:: eval_grid => eval_SphericalHarmonic_2Dgrid_GS, &
                                   eval_SphericalHarmonic_3Dgrid_GS
     procedure,private,pass:: calc_SphericalHarmonic_2Damps_GS
     procedure,private,pass:: calc_SphericalHarmonic_3Damps_GS
     procedure,private,pass:: eval_SphericalHarmonic_2Dgrid_GS
     procedure,private,pass:: eval_SphericalHarmonic_3Dgrid_GS
  end type SphericalHarmonic_GS_t

  type ZonalMean_t
     private
     integer             :: nbas
     integer             :: nring
     real(r8),allocatable:: area (:,:)
     real(r8),allocatable:: basis(:,:,:)
     real(r8),allocatable:: map  (:,:)
    contains
     procedure,pass:: init      => init_ZonalMean
     procedure,pass:: final     => final_ZonalMean
     generic,public:: calc_amps => calc_ZonalMean_2Damps, &
                                   calc_ZonalMean_3Damps
     generic,public:: eval_grid => eval_ZonalMean_2Dgrid, &
                                   eval_ZonalMean_3Dgrid
     procedure,private,pass:: calc_ZonalMean_2Damps
     procedure,private,pass:: calc_ZonalMean_3Damps
     procedure,private,pass:: eval_ZonalMean_2Dgrid
     procedure,private,pass:: eval_ZonalMean_3Dgrid
 end type ZonalMean_t

  type ZonalProfile_t
     private
     integer             :: nlat
     integer             :: nbas
     real(r8),allocatable:: area (:)
     real(r8),allocatable:: basis(:,:)
     real(r8),allocatable:: map  (:,:)
    contains
     procedure,pass:: init      => init_ZonalProfile
     procedure,pass:: final     => final_ZonalProfile
     generic,public:: calc_amps => calc_ZonalProfile_1Damps, &
                                   calc_ZonalProfile_2Damps
     generic,public:: eval_grid => eval_ZonalProfile_1Dgrid, &
                                   eval_ZonalProfile_2Dgrid
     procedure,private,pass:: calc_ZonalProfile_1Damps
     procedure,private,pass:: calc_ZonalProfile_2Damps
     procedure,private,pass:: eval_ZonalProfile_1Dgrid
     procedure,private,pass:: eval_ZonalProfile_2Dgrid
  end type ZonalProfile_t

  type ZonalAverage_t
     private
     logical             :: LINEARWGTS = .true.
     integer             :: nlat
     real(r8),allocatable:: area   (:)
     real(r8),allocatable:: a_norm (:)
     real(r8),allocatable:: area_g (:,:)
     integer ,allocatable:: idx_map(:,:)
     integer ,allocatable:: idx_wgt(:,:)
    contains
     procedure,pass:: init   => init_ZonalAverage
     procedure,pass:: final  => final_ZonalAverage
     generic,public:: binAvg => calc_ZonalAverage_2DbinAvg, &
                                calc_ZonalAverage_3DbinAvg
     generic,public:: setVal => set_ZonalAverage_2Dgrid, &
                                set_ZonalAverage_3Dgrid
     procedure,private,pass:: calc_ZonalAverage_2DbinAvg
     procedure,private,pass:: calc_ZonalAverage_3DbinAvg
     procedure,private,pass:: set_ZonalAverage_2Dgrid
     procedure,private,pass:: set_ZonalAverage_3Dgrid
  end type ZonalAverage_t

  real(r8), parameter :: halfPI = 0.5_r8*pi
  real(r8), parameter :: twoPI  = 2.0_r8*pi
  real(r8), parameter :: fourPI = 4.0_r8*pi
  real(r8), parameter :: qrtrPI = 0.25_r8*pi
  real(r8), parameter :: invSqrt4pi = 1._r8/sqrt(fourPI)

contains
    !=======================================================================
    subroutine init_SphericalHarmonic_COV(this,I_nmax,O_nbas,SAMPLE_NRING)
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
      class(SphericalHarmonic_COV_t) :: this
      integer ,intent(in ):: I_nmax
      integer ,intent(out):: O_nbas
      integer ,intent(in ),optional:: SAMPLE_NRING
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Clons(:,:)
      real(r8),allocatable:: Clats(:,:)
      real(r8),allocatable:: Bcoef(:)
      real(r8),allocatable:: Csum (:,:)
      real(r8),allocatable:: Cvec (:)
      real(r8),allocatable:: Bsum (:,:)
      real(r8),allocatable:: Bnorm(:)
      real(r8),allocatable:: Bcov (:,:)
      real(r8):: area(pcols),rlat,rlon

      real(r8),allocatable:: Nr    (:)
      real(r8),allocatable:: Nr_sum(:)
      real(r8),allocatable:: Rs    (:)
      real(r8),allocatable:: As    (:)
      real(r8):: Rc,Dth,Ws
      real(r8):: Slat0,Slon0,X0,Y0,Z0  
      real(r8):: Slon, Slat,Xs,Ys,Zs
      real(r8):: Rot(3,3),Rse,Bavg,Pnm
      integer :: kk,ns,Nsamp

      integer :: mm,nn,n2,nb,lchnk,ncols,cc
      integer :: cnum,Cvec_len

      integer :: nlcols, count, astat
      character(len=*), parameter :: subname = 'init_SphericalHarmonic_COV'

      if (I_nmax<1) then
         call endrun('SphericalHarmonic_COV%init: ERROR I_nmax must be greater than 0')
      end if

      ! Allocate space
      !-----------------
      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)
      if(allocated(this%map  )) deallocate(this%map)
      if(allocated(this%bindx)) deallocate(this%bindx)

      O_nbas = (I_nmax+1)**2
      this%nmax = I_nmax
      this%nbas = O_nbas
      allocate(this%area (pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'this%area')
      allocate(this%basis(pcols,begchunk:endchunk,O_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%basis')
      allocate(this%map  (O_nbas,O_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%map')
      allocate(this%bindx(0:O_nbas,-O_nbas:O_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%bindx')
      this%area (:,:)   = 0._r8
      this%basis(:,:,:) = 0._r8
      this%map  (:,:)   = 0._r8
      this%bindx(:,:)   = -1


      Cvec_len = 0
      do nn= 1,this%nbas
         do n2=nn,this%nbas
            Cvec_len = Cvec_len + 1
         end do
      end do

      nlcols = get_nlcols_p()

      allocate(Clons(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'Clons')
      allocate(Clats(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'Clats')
      allocate(Bcoef(O_nbas/2+1), stat=astat)
      call handle_allocate_error(astat, subname, 'Bcoef')
      allocate(Csum (nlcols, Cvec_len), stat=astat)
      call handle_allocate_error(astat, subname, 'Csum')
      allocate(Cvec (Cvec_len), stat=astat)
      call handle_allocate_error(astat, subname, 'Cvec')
      allocate(Bsum (nlcols, O_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bsum')
      allocate(Bnorm(O_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bnorm')
      allocate(Bcov (O_nbas,O_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bcov')

      Bsum(:,:) = 0._r8
      Csum(:,:) = 0._r8

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
            Clats    (cc,lchnk) = rlat
         end do
      end do

      ! Initialize Samping Grid Size
      !-----------------------------
      if(present(SAMPLE_NRING)) then
        this%nring = SAMPLE_NRING
      else
        this%nring = 1
      endif

      if(this%nring < 1) then
        call endrun('SphericalHarmonic_COV%init: ERROR nring be 1 or more')
      endif

      ! Set the number of gridpoints in each concentric ring
      !-----------------------------------------------------
      allocate(Nr    (this%nring))
      allocate(Nr_sum(this%nring))
      Nr    (1) = 1
      Nr_sum(1) = 1
      do kk=2,this%nring
        Nr    (kk) = 8*(kk-1)
        Nr_sum(kk) = Nr_sum(kk-1) + Nr(kk)
      end do
      Nsamp = Nr_sum(this%nring)

      ! Init the equal area grid with the center point, 
      ! then add the polar coordinate gridpoints for 
      ! each ring for a reference domain radius=1.
      !-------------------------------------------------
      allocate(Rs(Nsamp))
      allocate(As(Nsamp))
      Rs(1) = 0._r8
      As(1) = 0._r8
      ns = 1
      do kk =2,this%nring
        Rc = (sqrt(real(Nr_sum(kk)))+sqrt(real(Nr_sum(kk-1))))/(2._r8*sqrt(real(Nsamp)))
        Dth = twoPI/real(Nr(kk))
        do nn=1,Nr(kk)
          ns = ns + 1
          Rs(ns) = Rc
          As(ns) = Dth*(nn-1)
        end do
      end do
      Ws = 1._r8/real(Nsamp)

      if(ns.ne.Nsamp) then
        call endrun('SphericalHarmonic_COV%init: ERROR mismatch in sample point number')
      endif

      ! Samping INFO/DIAG Output
      !--------------------------
      IF((.TRUE.).and.(Nsamp>1)) THEN
        write(iulog,*) ' '
        write(iulog,*) 'SphericalHarmonic_COV%init: SAMPLE GRID: Nsamp=',Nsamp
        do ns=1,Nsamp
          write(iulog,*) ' ns=',ns,' Grid: Rs=',Rs(ns),' As=',As(ns),As(ns)*360./twoPI
        end do
        write(iulog,*) ' '
      ENDIF

      ! Add first basis for the mean values.
      !------------------------------------------
      this%nbas       = 1
      this%bindx(0,0) = 1
      this%basis(:,begchunk:endchunk,1) = invSqrt4pi

      ! Loop over the remaining meridional modes
      !------------------------------------------
      do nb=2,(this%nmax+1)
         nn = nb-1

         ! Loop over zonal modes for each nn
         !-----------------------------------
         do mm=-nn,nn
           this%nbas         = this%nbas + 1
           this%bindx(nn,mm) = this%nbas

           ! Generate coefs for the basis
           !------------------------------
           call sh_gen_basis_coefs(nn,iabs(mm),Bcoef)

           if(this%nring.eq.1) then
             ! Create basis for the coefs at each ncol gridpoint
             !    Lats shifted by PI/2 for SH's SP origin
             !---------------------------------------------------
             do lchnk=begchunk,endchunk
               ncols = get_ncols_p(lchnk)
               do cc = 1,ncols
                 Slat = Clats(cc,lchnk) + halfPI
                 call sh_create_basis(nn,iabs(mm),Slat,Bcoef,this%basis(cc,lchnk,this%nbas))
                 if(mm < 0) then
                   this%basis(cc,lchnk,this%nbas) = this%basis(cc,lchnk,this%nbas)*sin(mm*Clons(cc,lchnk))
                 elseif(mm > 0) then
                   this%basis(cc,lchnk,this%nbas) = this%basis(cc,lchnk,this%nbas)*cos(mm*Clons(cc,lchnk))
                 endif
               end do
             end do
           else
             ! Optionally Compute local area averge of 
             ! basis values for each gridpoiont
             !---------------------------------------
             do lchnk=begchunk,endchunk
               ncols = get_ncols_p(lchnk)
               do cc = 1,ncols

                 ! Set values to rotate/scale the sample 
                 ! grid to the current SE gridpoint
                 !-----------------------------------------
                 Rse      = sqrt(this%area(cc,lchnk)/pi)
                 Rot(1,1) = -sin(Clons(cc,lchnk))
                 Rot(2,1) = -cos(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
                 Rot(3,1) =  cos(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
                 Rot(1,2) =  cos(Clons(cc,lchnk))
                 Rot(2,2) = -sin(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
                 Rot(3,2) =  sin(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
                 Rot(1,3) =  0._r8
                 Rot(2,3) =  cos(Clats(cc,lchnk))
                 Rot(3,3) =  sin(Clats(cc,lchnk))

                 ! For each SE grid point, we rotate/scale the reference
                 ! grid for a local sampling domain with the elements area.
                 ! The computed basis values are then averaged over this domain
                 !------------------------------------------------------------
                 Bavg = 0._r8
                 do ns=1,Nsamp
                   ! scale polar coordinates to lat/lon circle centered at NP
                   !---------------------------------------------------------
                   Slat0 = halfPI - Rse*Rs(ns)
                   Slon0 = As(ns)
   
                   ! Compute cartesian coordinates
                   !-------------------------------
                   X0 = cos(Slon0)*cos(Slat0)
                   Y0 = sin(Slon0)*cos(Slat0)
                   Z0 = sin(Slat0)
   
                   ! Apply the Rotation from the NP
                   !--------------------------------
                   Xs = X0*Rot(1,1) + Y0*Rot(2,1) + Z0*Rot(3,1)
                   Ys = X0*Rot(1,2) + Y0*Rot(2,2) + Z0*Rot(3,2)
                   Zs = X0*Rot(1,3) + Y0*Rot(2,3) + Z0*Rot(3,3)
   
                   ! Compute resulting lat/lon gridpoint
                   !-------------------------------------
                   Slon = atan2(Ys,Xs)
                   Slat = halfPI + asin(Zs)
   
                   ! Compute basis value and add the result to the average
                   !    Slat shifted by PI/2 for SH's SP origin
                   !-------------------------------------------------------
                   Slat = Slat + halfPI
                   call sh_create_basis(nn,iabs(mm),Slat,Bcoef,Pnm)
                   if(mm < 0) then
                     Bavg = Bavg + Ws*Pnm*sin(mm*Slon)
                   elseif(mm > 0) then
                     Bavg = Bavg + Ws*Pnm*cos(mm*Slon)
                   endif
                 end do
                 this%basis(cc,lchnk,this%nbas) = Bavg
               end do
             end do
           endif

         end do ! mm=-nn,nn
      end do ! nn=2,this%nbas
      O_nbas = this%nbas

      ! Numerically normalize the basis funnctions
      !--------------------------------------------------------------
      do nn=1,this%nbas
         count = 0
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols
               count=count+1
               Bsum(count,nn) = this%basis(cc,lchnk,nn)*this%basis(cc,lchnk,nn)*this%area(cc,lchnk)
            end do
         end do
      end do ! nn=1,this%nbas

      call shr_reprosum_calc(Bsum, Bnorm, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

      do nn=1,this%nbas
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            this%basis(:ncols,lchnk,nn) = this%basis(:ncols,lchnk,nn)/sqrt(Bnorm(nn))
         end do
      end do ! nn=1,this%nbas

      ! Compute covariance matrix for basis functions
      ! (Yes, they are theoretically orthonormal, but lets make sure)
      !---------------------------------------------------------------
      cnum = 0
      do nn= 1,this%nbas
         do n2=nn,this%nbas
            cnum = cnum + 1
            count = 0
            do lchnk=begchunk,endchunk
               ncols = get_ncols_p(lchnk)
               do cc = 1,ncols
                  count=count+1
                  Csum(count,cnum) = this%basis(cc,lchnk,nn)*this%basis(cc,lchnk,n2)*this%area(cc,lchnk)
               end do
            end do

         end do
      end do

      call shr_reprosum_calc(Csum, Cvec, count, nlcols, Cvec_len, gbl_count=ngcols_p, commid=mpicom)

      cnum = 0
      do nn= 1,this%nbas
         do n2=nn,this%nbas
            cnum = cnum + 1
            Bcov(nn,n2) = Cvec(cnum)
            Bcov(n2,nn) = Cvec(cnum)
         end do
      end do

      ! Invert to get the basis amplitude map
      !--------------------------------------
      call Invert_Matrix(Bcov,this%nbas,this%map)

      ! End Routine
      !------------
      deallocate(Nr)
      deallocate(Nr_sum)
      deallocate(Rs)
      deallocate(As)
      deallocate(Clons)
      deallocate(Clats)
      deallocate(Bcoef)
      deallocate(Csum )
      deallocate(Cvec )
      deallocate(Bsum )
      deallocate(Bnorm)
      deallocate(Bcov )
    end subroutine init_SphericalHarmonic_COV
    !=======================================================================


    !=======================================================================
    subroutine final_SphericalHarmonic_COV(this)
      class(SphericalHarmonic_COV_t) :: this

      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)
      if(allocated(this%map  )) deallocate(this%map)

    end subroutine final_SphericalHarmonic_COV
    !=======================================================================


    !=======================================================================
    subroutine calc_SphericalHarmonic_2Damps_COV(this,I_Gdata,O_Bamp)
      !
      ! calc_SphericalHarmonic_2Damps: Given 2D data values for the ncol gridpoints,
      !                        compute the spherical harmonic basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_COV_t) :: this
      real(r8),intent(in ) :: I_Gdata(pcols,begchunk:endchunk)
      real(r8),intent(out) :: O_Bamp(:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable :: Csum(:,:)
      real(r8),allocatable :: Gcov(:)
      integer :: nn,n2,ncols,lchnk,cc
      integer :: nlcols, count, astat

      character(len=*), parameter :: subname = 'calc_SphericalHarmonic_2Damps_COV'

      nlcols = get_nlcols_p()

      allocate(Gcov(this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Gcov')
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

      call shr_reprosum_calc(Csum, Gcov, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

      ! Multiply by map to get the amplitudes
      !-------------------------------------------
      do nn=1,this%nbas
        O_Bamp(nn) = 0._r8
        do n2=1,this%nbas
          O_Bamp(nn) = O_Bamp(nn) + this%map(n2,nn)*Gcov(n2)
        end do
      end do

      ! End Routine
      !------------
      deallocate(Csum)
      deallocate(Gcov)

    end subroutine calc_SphericalHarmonic_2Damps_COV
    !=======================================================================


    !=======================================================================
    subroutine calc_SphericalHarmonic_3Damps_COV(this,I_Gdata,O_Bamp)
      !
      ! calc_SphericalHarmonic_3Damps: Given 3D data values for the ncol,nlev gridpoints,
      !                        compute the zonal mean basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_COV_t) :: this
      real(r8),intent(in ):: I_Gdata(:,:,begchunk:)
      real(r8),intent(out):: O_Bamp (:,:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Csum (:,:)
      real(r8),allocatable:: Gcov (:)
      integer:: nn,n2,ncols,lchnk,cc
      integer:: Nsum,ns,ll
      integer :: nlcols, count, astat

      integer :: nlev
      character(len=*), parameter :: subname = 'calc_SphericalHarmonic_3Damps_COV'

      nlev = size(I_Gdata,dim=2)

      nlcols = get_nlcols_p()
      allocate(Gcov(this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Gcov')
      allocate(Csum(nlcols, this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Csum')

      Csum(:,:) = 0._r8
      O_Bamp(:,:) = 0._r8

      ! Compute Covariance with input data and basis functions
      !--------------------------------------------------------
      do ll= 1,nlev

         Csum(:,:) = 0._r8
         Gcov(:) = 0._r8

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

         call shr_reprosum_calc(Csum, Gcov, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

         ! Multiply by map to get the amplitudes
         !-------------------------------------------
         do nn=1,this%nbas
            O_Bamp(nn,ll) = 0._r8
            do n2=1,this%nbas
               O_Bamp(nn,ll) = O_Bamp(nn,ll) + this%map(n2,nn)*Gcov(n2)
            end do
         end do

      end do

      ! End Routine
      !------------
      deallocate(Csum)
      deallocate(Gcov)

    end subroutine calc_SphericalHarmonic_3Damps_COV
    !=======================================================================


    !=======================================================================
    subroutine eval_SphericalHarmonic_2Dgrid_COV(this,I_Bamp,O_Gdata)
      !
      ! eval_SphericalHarmonic_2Dgrid: Given the zonal mean basis amplitudes,
      !                        compute 2D data values for the ncol gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_COV_t) :: this
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

    end subroutine eval_SphericalHarmonic_2Dgrid_COV
    !=======================================================================


    !=======================================================================
    subroutine eval_SphericalHarmonic_3Dgrid_COV(this,I_Bamp,O_Gdata)
      !
      ! eval_SphericalHarmonic_3Dgrid: Given the zonal mean basis amplitudes,
      !                      compute 3D data values for the ncol,nlev gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_COV_t) :: this
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

    end subroutine eval_SphericalHarmonic_3Dgrid_COV
    !=======================================================================


    !=======================================================================
    subroutine init_SphericalHarmonic_GS(this,I_nmax,O_nbas,SAMPLE_NRING)
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
      class(SphericalHarmonic_GS_t) :: this
      integer ,intent(in ):: I_nmax
      integer ,intent(out):: O_nbas
      integer ,intent(in ),optional:: SAMPLE_NRING
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

      real(r8),allocatable:: Nr    (:)
      real(r8),allocatable:: Nr_sum(:)
      real(r8),allocatable:: Rs    (:)
      real(r8),allocatable:: As    (:)
      real(r8):: Rc,Dth,Ws
      real(r8):: Slat0,Slon0,X0,Y0,Z0  
      real(r8):: Slon, Slat,Xs,Ys,Zs
      real(r8):: Rot(3,3),Rse,Bavg,Bavg_e,Bavg_o
      integer :: kk,ns,Nsamp

      integer :: nbas_e,nbas_o
      integer :: mm,nn,n2,nb,lchnk,ncols,cc

      integer :: nlcols, count, astat
      character(len=*), parameter :: subname = 'init_SphericalHarmonic_GS'

      if (I_nmax<1) then
         call endrun('SphericalHarmonic_GS%init: ERROR I_nmax must be greater than 0')
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
          Clats    (cc,lchnk) = rlat
        end do
      end do

      ! Initialize Samping Grid Size
      !-----------------------------
      if(present(SAMPLE_NRING)) then
        this%nring = SAMPLE_NRING
      else
        this%nring = 1
      endif

      if(this%nring < 1) then
        call endrun('SphericalHarmonic_GS%init: ERROR nring be 1 or more')
      endif

      ! Set the number of gridpoints in each concentric ring
      !-----------------------------------------------------
      allocate(Nr    (this%nring))
      allocate(Nr_sum(this%nring))
      Nr    (1) = 1
      Nr_sum(1) = 1
      do kk=2,this%nring
        Nr    (kk) = 8*(kk-1)
        Nr_sum(kk) = Nr_sum(kk-1) + Nr(kk)
      end do
      Nsamp = Nr_sum(this%nring)

      ! Init the equal area grid with the center point, 
      ! then add the polar coordinate gridpoints for 
      ! each ring for a reference domain radius=1.
      !-------------------------------------------------
      allocate(Rs(Nsamp))
      allocate(As(Nsamp))
      Rs(1) = 0._r8
      As(1) = 0._r8
      ns = 1
      do kk =2,this%nring
        Rc = (sqrt(real(Nr_sum(kk)))+sqrt(real(Nr_sum(kk-1))))/(2._r8*sqrt(real(Nsamp)))
        Dth = twoPI/real(Nr(kk))
        do nn=1,Nr(kk)
          ns = ns + 1
          Rs(ns) = Rc
          As(ns) = Dth*(nn-1)
        end do
      end do
      Ws = 1._r8/real(Nsamp)

      if(ns.ne.Nsamp) then
        call endrun('SphericalHarmonic_GS%init: ERROR mismatch in sample point number')
      endif

      ! Samping INFO/DIAG Output
      !--------------------------
      IF((.TRUE.).and.(Nsamp>1)) THEN
        write(iulog,*) ' '
        write(iulog,*) 'SphericalHarmonic_GS%init: SAMPLE GRID: Nsamp=',Nsamp
        do ns=1,Nsamp
          write(iulog,*) ' ns=',ns,' Grid: Rs=',Rs(ns),' As=',As(ns),As(ns)*360./twoPI
        end do
        write(iulog,*) ' '
      ENDIF

      ! Add first basis for the mean values.
      !------------------------------------------
      this%nbas = 1
      this%basis(:,begchunk:endchunk,1) = invSqrt4pi

      ! Loop over the remaining meridional modes, 
      ! The ordering of the basis functions here is 
      ! set to accomodate Gram-Schmidt
      !------------------------------------------------
      do nb=2,(this%nmax+1)
        nn = nb-1

        ! Add the m=0 mode first
        !------------------------
        this%nbas = this%nbas + 1
        call sh_gen_basis_coefs(nn,0,Bcoef)

        if(this%nring.eq.1) then
          do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols
              !  Lats shifted by PI/2 for SH's SP origin
              !--------------------------------------------
              Slat = Clats(cc,lchnk) + halfPI
              call sh_create_basis(nn,0,Slat,Bcoef,this%basis(cc,lchnk,this%nbas))
            end do
          end do
        else
          ! Optionally Compute local area averge of 
          ! basis values for each gridpoiont
          !---------------------------------------
          do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols

              ! Set values to rotate/scale the sample 
              ! grid to the current SE gridpoint
              !-----------------------------------------
              Rse      = sqrt(this%area(cc,lchnk)/pi)
              Rot(1,1) = -sin(Clons(cc,lchnk))
              Rot(2,1) = -cos(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
              Rot(3,1) =  cos(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
              Rot(1,2) =  cos(Clons(cc,lchnk))
              Rot(2,2) = -sin(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
              Rot(3,2) =  sin(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
              Rot(1,3) =  0._r8
              Rot(2,3) =  cos(Clats(cc,lchnk))
              Rot(3,3) =  sin(Clats(cc,lchnk))

              ! For each SE grid point, we rotate/scale the reference
              ! grid for a local sampling domain with the elements area.
              ! The computed basis values are then averaged over this domain
              !------------------------------------------------------------
              Bavg = 0._r8
              do ns=1,Nsamp
                ! scale polar coordinates to lat/lon circle centered at NP
                !---------------------------------------------------------
                Slat0 = halfPI - Rse*Rs(ns)
                Slon0 = As(ns)
   
                ! Compute cartesian coordinates
                !-------------------------------
                X0 = cos(Slon0)*cos(Slat0)
                Y0 = sin(Slon0)*cos(Slat0)
                Z0 = sin(Slat0)
   
                ! Apply the Rotation from the NP
                !--------------------------------
                Xs = X0*Rot(1,1) + Y0*Rot(2,1) + Z0*Rot(3,1)
                Ys = X0*Rot(1,2) + Y0*Rot(2,2) + Z0*Rot(3,2)
                Zs = X0*Rot(1,3) + Y0*Rot(2,3) + Z0*Rot(3,3)
   
                ! Compute resulting lat/lon gridpoint
                !-------------------------------------
                Slon = atan2(Ys,Xs)
                Slat = halfPI + asin(Zs)
   
                ! Compute basis value and add the result to the average
                !    Slat shifted by PI/2 for SH's SP origin
                !-------------------------------------------------------
                Slat = Slat + halfPI
                call sh_create_basis(nn,0,Slat,Bcoef,Pnm)
                Bavg = Bavg + Ws*Pnm
              end do
              this%basis(cc,lchnk,this%nbas) = Bavg
            end do
          end do
        endif

        ! Now loop over zonal modes mm=1,nn
        ! and add even/odd basis functions
        !----------------------------------
        do mm=1,nn
          nbas_e = this%nbas + 1
          nbas_o = this%nbas + 2
          this%nbas = this%nbas + 2
          call sh_gen_basis_coefs(nn,mm,Bcoef)

          if(this%nring.eq.1) then
            do lchnk=begchunk,endchunk
              ncols = get_ncols_p(lchnk)
              do cc = 1,ncols
                !  Lats shifted by PI/2 for SH's SP origin
                !--------------------------------------------
                Slat = Clats(cc,lchnk) + halfPI
                call sh_create_basis(nn,mm,Slat,Bcoef,Pnm)
                this%basis(cc,lchnk,nbas_e) = Pnm*cos(mm*Clons(cc,lchnk))
                this%basis(cc,lchnk,nbas_o) = Pnm*sin(mm*Clons(cc,lchnk))
              end do
            end do
          else
            ! Optionally Compute local area averge of 
            ! basis values for each gridpoiont
            !---------------------------------------
            do lchnk=begchunk,endchunk
              ncols = get_ncols_p(lchnk)
              do cc = 1,ncols

                ! Set values to rotate/scale the sample 
                ! grid to the current SE gridpoint
                !-----------------------------------------
                Rse      = sqrt(this%area(cc,lchnk)/pi)
                Rot(1,1) = -sin(Clons(cc,lchnk))
                Rot(2,1) = -cos(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
                Rot(3,1) =  cos(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
                Rot(1,2) =  cos(Clons(cc,lchnk))
                Rot(2,2) = -sin(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
                Rot(3,2) =  sin(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
                Rot(1,3) =  0._r8
                Rot(2,3) =  cos(Clats(cc,lchnk))
                Rot(3,3) =  sin(Clats(cc,lchnk))

                ! For each SE grid point, we rotate/scale the reference
                ! grid for a local sampling domain with the elements area.
                ! The computed basis values are then averaged over this domain
                !------------------------------------------------------------
                Bavg_e = 0._r8
                Bavg_o = 0._r8
                do ns=1,Nsamp
                  ! scale polar coordinates to lat/lon circle centered at NP
                  !---------------------------------------------------------
                  Slat0 = halfPI - Rse*Rs(ns)
                  Slon0 = As(ns)
   
                  ! Compute cartesian coordinates
                  !-------------------------------
                  X0 = cos(Slon0)*cos(Slat0)
                  Y0 = sin(Slon0)*cos(Slat0)
                  Z0 = sin(Slat0)
   
                  ! Apply the Rotation from the NP
                  !--------------------------------
                  Xs = X0*Rot(1,1) + Y0*Rot(2,1) + Z0*Rot(3,1)
                  Ys = X0*Rot(1,2) + Y0*Rot(2,2) + Z0*Rot(3,2)
                  Zs = X0*Rot(1,3) + Y0*Rot(2,3) + Z0*Rot(3,3)
   
                  ! Compute resulting lat/lon gridpoint
                  !-------------------------------------
                  Slon = atan2(Ys,Xs)
                  Slat = halfPI + asin(Zs)
   
                  ! Compute basis value and add the result to the average
                  !    Slat shifted by PI/2 for SH's SP origin
                  !-------------------------------------------------------
                  Slat = Slat + halfPI
                  call sh_create_basis(nn,mm,Slat,Bcoef,Pnm)
                  Bavg_o = Bavg_o + Ws*Pnm*sin(mm*Slon)
                  Bavg_e = Bavg_e + Ws*Pnm*cos(mm*Slon)
                end do
                this%basis(cc,lchnk,nbas_e) = Bavg_e
                this%basis(cc,lchnk,nbas_o) = Bavg_o
              end do
            end do
          endif

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
      deallocate(Nr)
      deallocate(Nr_sum)
      deallocate(Rs)
      deallocate(As)
      deallocate(Clons)
      deallocate(Clats)
      deallocate(Bcoef)
      deallocate(Bsum )
      deallocate(Bamp )
      deallocate(Bnorm)
    end subroutine init_SphericalHarmonic_GS
    !=======================================================================


    !=======================================================================
    subroutine final_SphericalHarmonic_GS(this)
      class(SphericalHarmonic_GS_t) :: this

      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)

    end subroutine final_SphericalHarmonic_GS
    !=======================================================================


    !=======================================================================
    subroutine calc_SphericalHarmonic_2Damps_GS(this,I_Gdata,O_Bamp)
      !
      ! calc_SphericalHarmonic_2Damps: Given 2D data values for the ncol gridpoints,
      !                        compute the spherical harmonic basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_GS_t) :: this
      real(r8),intent(in ) :: I_Gdata(pcols,begchunk:endchunk)
      real(r8),intent(out) :: O_Bamp(:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable :: Csum(:,:)
      real(r8),allocatable :: Bamp(:)
      integer :: nn,n2,ncols,lchnk,cc
      integer :: nlcols, count, astat

      character(len=*), parameter :: subname = 'calc_SphericalHarmonic_2Damps_GS'

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

    end subroutine calc_SphericalHarmonic_2Damps_GS
    !=======================================================================


    !=======================================================================
    subroutine calc_SphericalHarmonic_3Damps_GS(this,I_Gdata,O_Bamp)
      !
      ! calc_SphericalHarmonic_3Damps: Given 3D data values for the ncol,nlev gridpoints,
      !                        compute the zonal mean basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_GS_t) :: this
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
      character(len=*), parameter :: subname = 'calc_SphericalHarmonic_3Damps_GS'

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

    end subroutine calc_SphericalHarmonic_3Damps_GS
    !=======================================================================


    !=======================================================================
    subroutine eval_SphericalHarmonic_2Dgrid_GS(this,I_Bamp,O_Gdata)
      !
      ! eval_SphericalHarmonic_2Dgrid: Given the zonal mean basis amplitudes,
      !                        compute 2D data values for the ncol gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_GS_t) :: this
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

    end subroutine eval_SphericalHarmonic_2Dgrid_GS
    !=======================================================================


    !=======================================================================
    subroutine eval_SphericalHarmonic_3Dgrid_GS(this,I_Bamp,O_Gdata)
      !
      ! eval_SphericalHarmonic_3Dgrid: Given the zonal mean basis amplitudes,
      !                      compute 3D data values for the ncol,nlev gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(SphericalHarmonic_GS_t) :: this
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

    end subroutine eval_SphericalHarmonic_3Dgrid_GS
    !=======================================================================


    !=======================================================================
    subroutine init_ZonalMean(this,I_nbas,SAMPLE_NRING)
      !
      ! init_ZonalMean: Initialize the ZonalMean data structures for the
      !                 physics grid. It is assumed that the domain
      !                 of these gridpoints spans the surface of the sphere.
      !                 The representation of basis functions is
      !                 normalized w.r.t integration over the sphere.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalMean_t) :: this
      integer ,intent(in):: I_nbas
      integer ,intent(in),optional:: SAMPLE_NRING
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Clons(:,:)
      real(r8),allocatable:: Clats(:,:)
      real(r8),allocatable:: Bcoef(:)
      real(r8),allocatable:: Csum (:,:)
      real(r8),allocatable:: Cvec (:)
      real(r8),allocatable:: Bsum (:,:)
      real(r8),allocatable:: Bnorm(:)
      real(r8),allocatable:: Bcov (:,:)
      real(r8):: area(pcols),rlat,rlon

      real(r8),allocatable:: Nr    (:)
      real(r8),allocatable:: Nr_sum(:)
      real(r8),allocatable:: Rs    (:)
      real(r8),allocatable:: As    (:)
      real(r8):: Rc,Dth,Ws
      real(r8):: Slat0,Slon0,X0,Y0,Z0  
      real(r8):: Slon, Slat,Xs,Ys,Zs
      real(r8):: Rot(3,3),Rse,Bavg,Pnm
      integer :: kk,ns,mm,Nsamp

      integer :: nn,n2,nb,lchnk,ncols,cc
      integer :: cnum,Cvec_len

      integer :: nlcols, count, astat
      character(len=*), parameter :: subname = 'init_ZonalMean'

      if (I_nbas<1) then
         call endrun('ZonalMean%init: ERROR I_nbas must be greater than 0')
      end if

      ! Allocate space
      !-----------------
      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)
      if(allocated(this%map  )) deallocate(this%map)

      this%nbas = I_nbas
      allocate(this%area (pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'this%area')
      allocate(this%basis(pcols,begchunk:endchunk,I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%basis')
      allocate(this%map  (I_nbas,I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%map')
      this%area (:,:)   = 0._r8
      this%basis(:,:,:) = 0._r8
      this%map  (:,:)   = 0._r8

      Cvec_len = 0
      do nn= 1,this%nbas
         do n2=nn,this%nbas
            Cvec_len = Cvec_len + 1
         end do
      end do

      nlcols = get_nlcols_p()

      allocate(Clons(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'Clons')
      allocate(Clats(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'Clats')
      allocate(Bcoef(I_nbas/2+1), stat=astat)
      call handle_allocate_error(astat, subname, 'Bcoef')
      allocate(Csum (nlcols, Cvec_len), stat=astat)
      call handle_allocate_error(astat, subname, 'Csum')
      allocate(Cvec (Cvec_len), stat=astat)
      call handle_allocate_error(astat, subname, 'Cvec')
      allocate(Bsum (nlcols, I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bsum')
      allocate(Bnorm(I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bnorm')
      allocate(Bcov (I_nbas,I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bcov')

      Bsum(:,:) = 0._r8
      Csum(:,:) = 0._r8

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
            Clats    (cc,lchnk) = rlat
            Clons    (cc,lchnk) = rlon
         end do
      end do

      ! Initialize Samping Grid Size
      !-----------------------------
      if(present(SAMPLE_NRING)) then
        this%nring = SAMPLE_NRING
      else
        this%nring = 1
      endif

      if(this%nring < 1) then
        call endrun('ZonalMean%init: ERROR nring be 1 or more')
      endif

      ! Set the number of gridpoints in each concentric ring
      !-----------------------------------------------------
      allocate(Nr    (this%nring))
      allocate(Nr_sum(this%nring))
      Nr    (1) = 1
      Nr_sum(1) = 1
      do kk=2,this%nring
        Nr    (kk) = 8*(kk-1)
        Nr_sum(kk) = Nr_sum(kk-1) + Nr(kk)
      end do
      Nsamp = Nr_sum(this%nring)

      ! Init the equal area grid with the center point, 
      ! then add the polar coordinate gridpoints for 
      ! each ring for a reference domain radius=1.
      !-------------------------------------------------
      allocate(Rs(Nsamp))
      allocate(As(Nsamp))
      Rs(1) = 0._r8
      As(1) = 0._r8
      ns = 1
      do kk =2,this%nring
        Rc = (sqrt(real(Nr_sum(kk)))+sqrt(real(Nr_sum(kk-1))))/(2._r8*sqrt(real(Nsamp)))
        Dth = twoPI/real(Nr(kk))
        do nn=1,Nr(kk)
          ns = ns + 1
          Rs(ns) = Rc
          As(ns) = Dth*(nn-1)
        end do
      end do
      Ws = 1._r8/real(Nsamp)

      if(ns.ne.Nsamp) then
        call endrun('ZonalMean%init: ERROR mismatch in sample point number')
      endif

      ! Samping INFO/DIAG Output
      !--------------------------
      IF((.TRUE.).and.(Nsamp>1)) THEN
        write(iulog,*) ' '
        write(iulog,*) 'ZonalMean%init: SAMPLE GRID: Nsamp=',Nsamp
        do ns=1,Nsamp
          write(iulog,*) ' ns=',ns,' Grid: Rs=',Rs(ns),' As=',As(ns),As(ns)*360./twoPI
        end do
        write(iulog,*) ' '
      ENDIF

      ! Add first basis for the mean values.
      !------------------------------------------
      this%basis(:,begchunk:endchunk,1) = invSqrt4pi

      ! Loop over the remaining basis functions
      !---------------------------------------
      do nn=2,this%nbas
         nb = nn-1

         ! Generate coefs for the basis
         !------------------------------
         call sh_gen_basis_coefs(nb,0,Bcoef)

         if(this%nring.eq.1) then
           ! Create basis for the coefs at each ncol gridpoint
           !    Lats shifted by PI/2 for SH's SP origin
           !---------------------------------------------------
           do lchnk=begchunk,endchunk
              ncols = get_ncols_p(lchnk)
              do cc = 1,ncols
                 Slat = Clats(cc,lchnk) + halfPI
                 call sh_create_basis(nb,0,Slat,Bcoef,this%basis(cc,lchnk,nn))
              end do
           end do
         else
           ! Optionally Compute local area averge of 
           ! basis values for each gridpoiont
           !---------------------------------------
           do lchnk=begchunk,endchunk
             ncols = get_ncols_p(lchnk)
             do cc = 1,ncols

               ! Set values to rotate/scale the sample 
               ! grid to the current SE gridpoint
               !-----------------------------------------
               Rse      = sqrt(this%area(cc,lchnk)/pi)
               Rot(1,1) = -sin(Clons(cc,lchnk))
               Rot(2,1) = -cos(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
               Rot(3,1) =  cos(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
               Rot(1,2) =  cos(Clons(cc,lchnk))
               Rot(2,2) = -sin(Clons(cc,lchnk))*sin(Clats(cc,lchnk))
               Rot(3,2) =  sin(Clons(cc,lchnk))*cos(Clats(cc,lchnk))
               Rot(1,3) =  0._r8
               Rot(2,3) =  cos(Clats(cc,lchnk))
               Rot(3,3) =  sin(Clats(cc,lchnk))

               ! For each SE grid point, we rotate/scale the reference
               ! grid for a local sampling domain with the elements area.
               ! The computed basis values are then averaged over this domain
               !------------------------------------------------------------
               Bavg = 0._r8
               do ns=1,Nsamp
                 ! scale polar coordinates to lat/lon circle centered at NP
                 !---------------------------------------------------------
                 Slat0 = halfPI - Rse*Rs(ns)
                 Slon0 = As(ns)
 
                 ! Compute cartesian coordinates
                 !-------------------------------
                 X0 = cos(Slon0)*cos(Slat0)
                 Y0 = sin(Slon0)*cos(Slat0)
                 Z0 = sin(Slat0)
 
                 ! Apply the Rotation from the NP
                 !--------------------------------
                 Xs = X0*Rot(1,1) + Y0*Rot(2,1) + Z0*Rot(3,1)
                 Ys = X0*Rot(1,2) + Y0*Rot(2,2) + Z0*Rot(3,2)
                 Zs = X0*Rot(1,3) + Y0*Rot(2,3) + Z0*Rot(3,3)
 
                 ! Compute resulting lat/lon gridpoint
                 !-------------------------------------
                 Slon = atan2(Ys,Xs)
                 Slat = halfPI + asin(Zs)
 
                 ! Compute basis value and add the result to the average
                 !    Slat shifted by PI/2 for SH's SP origin
                 !-------------------------------------------------------
                 Slat = Slat + halfPI
                 call sh_create_basis(nn,0,Slat,Bcoef,Pnm)
                 if(mm < 0) then
                   Bavg = Bavg + Ws*Pnm*sin(mm*Slon)
                 elseif(mm > 0) then
                   Bavg = Bavg + Ws*Pnm*cos(mm*Slon)
                 endif
               end do
               this%basis(cc,lchnk,this%nbas) = Bavg
             end do
           end do
         endif

      end do ! nn=2,this%nbas

      ! Numerically normalize the basis funnctions
      !--------------------------------------------------------------
      do nn=1,this%nbas
         count = 0
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols
               count=count+1
               Bsum(count,nn) = this%basis(cc,lchnk,nn)*this%basis(cc,lchnk,nn)*this%area(cc,lchnk)
            end do
         end do
      end do ! nn=1,this%nbas

      call shr_reprosum_calc(Bsum, Bnorm, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

      do nn=1,this%nbas
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            this%basis(:ncols,lchnk,nn) = this%basis(:ncols,lchnk,nn)/sqrt(Bnorm(nn))
         end do
      end do ! nn=1,this%nbas

      ! Compute covariance matrix for basis functions
      ! (Yes, they are theoretically orthonormal, but lets make sure)
      !---------------------------------------------------------------
      cnum = 0
      do nn= 1,this%nbas
         do n2=nn,this%nbas
            cnum = cnum + 1
            count = 0
            do lchnk=begchunk,endchunk
               ncols = get_ncols_p(lchnk)
               do cc = 1,ncols
                  count=count+1
                  Csum(count,cnum) = this%basis(cc,lchnk,nn)*this%basis(cc,lchnk,n2)*this%area(cc,lchnk)
               end do
            end do

         end do
      end do

      call shr_reprosum_calc(Csum, Cvec, count, nlcols, Cvec_len, gbl_count=ngcols_p, commid=mpicom)

      cnum = 0
      do nn= 1,this%nbas
         do n2=nn,this%nbas
            cnum = cnum + 1
            Bcov(nn,n2) = Cvec(cnum)
            Bcov(n2,nn) = Cvec(cnum)
         end do
      end do

      ! Invert to get the basis amplitude map
      !--------------------------------------
      call Invert_Matrix(Bcov,this%nbas,this%map)

      ! End Routine
      !------------
      deallocate(Nr)
      deallocate(Nr_sum)
      deallocate(Rs)
      deallocate(As)
      deallocate(Clats)
      deallocate(Bcoef)
      deallocate(Csum )
      deallocate(Cvec )
      deallocate(Bsum )
      deallocate(Bnorm)
      deallocate(Bcov )

    end subroutine init_ZonalMean
    !=======================================================================


    !=======================================================================
    subroutine final_ZonalMean(this)
      class(ZonalMean_t) :: this

      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)
      if(allocated(this%map  )) deallocate(this%map)

    end subroutine final_ZonalMean
    !=======================================================================


    !=======================================================================
    subroutine calc_ZonalMean_2Damps(this,I_Gdata,O_Bamp)
      !
      ! calc_ZonalMean_2Damps: Given 2D data values for the ncol gridpoints,
      !                        compute the zonal mean basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalMean_t) :: this
      real(r8),intent(in ) :: I_Gdata(pcols,begchunk:endchunk)
      real(r8),intent(out) :: O_Bamp(:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable :: Csum(:,:)
      real(r8),allocatable :: Gcov(:)
      integer :: nn,n2,ncols,lchnk,cc
      integer :: nlcols, count, astat

      character(len=*), parameter :: subname = 'calc_ZonalMean_2Damps'

      nlcols = get_nlcols_p()

      allocate(Gcov(this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Gcov')
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

      call shr_reprosum_calc(Csum, Gcov, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

      ! Multiply by map to get the amplitudes
      !-------------------------------------------
      do nn=1,this%nbas
        O_Bamp(nn) = 0._r8
        do n2=1,this%nbas
          O_Bamp(nn) = O_Bamp(nn) + this%map(n2,nn)*Gcov(n2)
        end do
      end do

      ! End Routine
      !------------
      deallocate(Csum)
      deallocate(Gcov)

    end subroutine calc_ZonalMean_2Damps
    !=======================================================================


    !=======================================================================
    subroutine calc_ZonalMean_3Damps(this,I_Gdata,O_Bamp)
      !
      ! calc_ZonalMean_3Damps: Given 3D data values for the ncol,nlev gridpoints,
      !                        compute the zonal mean basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalMean_t) :: this
      real(r8),intent(in ):: I_Gdata(:,:,begchunk:)
      real(r8),intent(out):: O_Bamp (:,:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Csum (:,:)
      real(r8),allocatable:: Gcov (:)
      integer:: nn,n2,ncols,lchnk,cc
      integer:: Nsum,ns,ll
      integer :: nlcols, count, astat

      integer :: nlev
      character(len=*), parameter :: subname = 'calc_ZonalMean_3Damps'

      nlev = size(I_Gdata,dim=2)

      nlcols = get_nlcols_p()
      allocate(Gcov(this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Gcov')
      allocate(Csum(nlcols, this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Csum')

      Csum(:,:) = 0._r8
      O_Bamp(:,:) = 0._r8

      ! Compute Covariance with input data and basis functions
      !--------------------------------------------------------
      do ll= 1,nlev

         Csum(:,:) = 0._r8
         Gcov(:) = 0._r8

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

         call shr_reprosum_calc(Csum, Gcov, count, nlcols, this%nbas, gbl_count=ngcols_p, commid=mpicom)

         ! Multiply by map to get the amplitudes
         !-------------------------------------------
         do nn=1,this%nbas
            O_Bamp(nn,ll) = 0._r8
            do n2=1,this%nbas
               O_Bamp(nn,ll) = O_Bamp(nn,ll) + this%map(n2,nn)*Gcov(n2)
            end do
         end do

      end do

      ! End Routine
      !------------
      deallocate(Csum)
      deallocate(Gcov)

    end subroutine calc_ZonalMean_3Damps
    !=======================================================================


    !=======================================================================
    subroutine eval_ZonalMean_2Dgrid(this,I_Bamp,O_Gdata)
      !
      ! eval_ZonalMean_2Dgrid: Given the zonal mean basis amplitudes,
      !                        compute 2D data values for the ncol gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalMean_t) :: this
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

    end subroutine eval_ZonalMean_2Dgrid
    !=======================================================================


    !=======================================================================
    subroutine eval_ZonalMean_3Dgrid(this,I_Bamp,O_Gdata)
      !
      ! eval_ZonalMean_3Dgrid: Given the zonal mean basis amplitudes,
      !                      compute 3D data values for the ncol,nlev gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalMean_t) :: this
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

    end subroutine eval_ZonalMean_3Dgrid
    !=======================================================================


    !=======================================================================
    subroutine init_ZonalProfile(this,IO_lats,IO_area,I_nlat,I_nbas,GEN_GAUSSLATS)
      !
      ! init_ZonalProfile: Initialize the ZonalProfile data structure for the
      !                    given nlat gridpoints. It is assumed that the domain
      !                    of these gridpoints of the profile span latitudes
      !                    from SP to NP.
      !                    The representation of basis functions functions is
      !                    normalized w.r.t integration over the sphere so that
      !                    when configured for tha same number of basis functions,
      !                    the calculated amplitudes are interchangable with
      !                    those for the ZonalMean_t class.
      !
      !                    The optional GEN_GAUSSLATS flag allows for the
      !                    generation of Gaussian latitudes. The generated grid
      !                    over-writes the values of IO_lats/IO_area passed by
      !                    the user.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalProfile_t) :: this
      real(r8)     ,intent(inout):: IO_lats(:)
      real(r8)     ,intent(inout):: IO_area(:)
      integer         ,intent(in):: I_nlat
      integer         ,intent(in):: I_nbas
      logical,optional,intent(in):: GEN_GAUSSLATS
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Clats(:)
      real(r8),allocatable:: Bcoef(:)
      real(r8),allocatable:: Bcov (:,:)
      real(r8):: Bnorm
      integer :: ii,nn,n2,nb,ierr, astat
      logical :: generate_lats

      character(len=*), parameter :: subname = 'init_ZonalProfile'

      generate_lats = .false.

      if (present(GEN_GAUSSLATS)) then
         generate_lats = GEN_GAUSSLATS
      end if

      ! Allocate space
      !-----------------
      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)
      if(allocated(this%map  )) deallocate(this%map)

      this%nlat = I_nlat
      this%nbas = I_nbas
      allocate(this%area (I_nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'this%area')
      allocate(this%basis(I_nlat,I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%basis')
      allocate(this%map  (I_nbas,I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'this%map')

      allocate(Clats(I_nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'Clats')
      allocate(Bcoef(I_nbas/2+1), stat=astat)
      call handle_allocate_error(astat, subname, 'Bcoef')
      allocate(Bcov (I_nbas,I_nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Bcov')

      ! Optionally create the Latitude Gridpoints
      ! and their associated area weights. Otherwise
      ! they need to be supplied by the user.
      !-----------------------------------------------
      if(generate_lats) then

        ! Create a Gaussian grid from SP to NP
        !--------------------------------------
        call sh_create_gaus_grid(I_nlat,Clats,IO_area,ierr)
        if (ierr/=0) then
           call endrun('init_ZonalProfile: Error creating Gaussian grid')
        end if

        ! Convert generated colatitudes SP->NP to Lats and convert
        ! to degrees and scale the area for global 2D integrals
        !-----------------------------------------------------------
        do nn=1,I_nlat
          IO_lats(nn) = (45._r8*Clats(nn)/qrtrPI) - 90._r8
          IO_area(nn) = IO_area(nn)*twoPI
        end do
      else
        ! Convert Latitudes to SP->NP colatitudes in radians
        !----------------------------------------------------
        do nn=1,I_nlat
          Clats(nn) = (IO_lats(nn) + 90._r8)*qrtrPI/45._r8
        end do
      endif

      ! Copy the area weights for each nlat
      ! gridpoint to the data structure
      !---------------------------------------
      this%area(1:I_nlat) = IO_area(1:I_nlat)

      ! Add first basis for the mean values.
      !------------------------------------------
      this%basis(:,1) = invSqrt4pi
      Bnorm = 0._r8
      do ii=1,I_nlat
        Bnorm = Bnorm + (this%basis(ii,1)*this%basis(ii,1)*this%area(ii))
      end do
      this%basis(:,1) = this%basis(:,1)/sqrt(Bnorm)

      ! Loop over the remaining basis functions
      !---------------------------------------
      do nn=2,I_nbas
        nb = nn-1

        ! Generate coefs for the basis
        !------------------------------
        call sh_gen_basis_coefs(nb,0,Bcoef)

        ! Create an un-normalized basis for the
        ! coefs at each nlat gridpoint
        !---------------------------------------
        do ii=1,I_nlat
          call sh_create_basis(nb,0,Clats(ii),Bcoef,this%basis(ii,nn))
        end do

        ! Numerically normalize the basis funnction
        !--------------------------------------------------------------
        Bnorm = 0._r8
        do ii=1,I_nlat
          Bnorm = Bnorm + (this%basis(ii,nn)*this%basis(ii,nn)*this%area(ii))
        end do
        this%basis(:,nn) = this%basis(:,nn)/sqrt(Bnorm)

      end do ! nn=1,I_nbas

      ! Compute covariance matrix for basis functions
      ! (Yes, they are theoretically orthonormal, but lets make sure)
      !--------------------------------------------------------------
      do nn=1,I_nbas
        do n2=1,I_nbas
          Bcov(nn,n2) = 0._r8
          do ii=1,I_nlat
            Bcov(nn,n2) = Bcov(nn,n2) + (this%basis(ii,nn)*this%basis(ii,n2)*this%area(ii))
          end do
        end do
      end do

      ! Invert to get the basis amplitude map
      !--------------------------------------
      call Invert_Matrix(Bcov,I_nbas,this%map)

      ! End Routine
      !------------
      deallocate(Clats)
      deallocate(Bcoef)
      deallocate(Bcov )

    end subroutine init_ZonalProfile
    !=======================================================================


    !=======================================================================
    subroutine final_ZonalProfile(this)
      class(ZonalProfile_t) :: this

      if(allocated(this%area )) deallocate(this%area)
      if(allocated(this%basis)) deallocate(this%basis)
      if(allocated(this%map  )) deallocate(this%map)

    end subroutine final_ZonalProfile
    !=======================================================================


    !=======================================================================
    subroutine calc_ZonalProfile_1Damps(this,I_Zdata,O_Bamp)
      !
      ! calc_ZonalProfile_1Damps: Given 1D data values for the nlat zonal
      !                           profiles gridpoints, compute the zonal
      !                           profile basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalProfile_t):: this
      real(r8),intent(in ):: I_Zdata(:)
      real(r8),intent(out):: O_Bamp (:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Gcov(:)
      integer:: ii,nn,n2, astat
      character(len=*), parameter :: subname = 'calc_ZonalProfile_1Damps'

      ! Compute Covariance with input data and basis functions
      !--------------------------------------------------------
      allocate(Gcov(this%nbas), stat=astat)
      call handle_allocate_error(astat, subname, 'Gcov')
      do nn=1,this%nbas
        Gcov(nn) = 0._r8
        do ii=1,this%nlat
          Gcov(nn) = Gcov(nn) + (I_Zdata(ii)*this%basis(ii,nn)*this%area(ii))
        end do
      end do

      ! Multiply by map to get the amplitudes
      !-------------------------------------------
      do nn=1,this%nbas
        O_Bamp(nn) = 0._r8
        do n2=1,this%nbas
          O_Bamp(nn) = O_Bamp(nn) + this%map(n2,nn)*Gcov(n2)
        end do
      end do

      deallocate(Gcov)

    end subroutine calc_ZonalProfile_1Damps
    !=======================================================================


    !=======================================================================
    subroutine calc_ZonalProfile_2Damps(this,I_Zdata,O_Bamp)
      !
      ! calc_ZonalProfile_2Damps: Given 2D data values for the nlat,nlev zonal
      !                           profiles gridpoints, compute the zonal
      !                           profile basis amplitudes.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalProfile_t):: this
      real(r8),intent(in ):: I_Zdata(:,:)
      real(r8),intent(out):: O_Bamp (:,:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Gcov(:,:)
      integer:: ii,nn,n2,ilev

      integer :: nlev, astat
      character(len=*), parameter :: subname = 'calc_ZonalProfile_2Damps'

      nlev = size(I_Zdata,dim=2)

      ! Compute Covariance with input data and basis functions
      !--------------------------------------------------------
      allocate(Gcov(this%nbas,nlev), stat=astat)
      call handle_allocate_error(astat, subname, 'Gcov')
      do ilev=1,nlev
         do nn=1,this%nbas
            Gcov(nn,ilev) = 0._r8
            do ii=1,this%nlat
               Gcov(nn,ilev) = Gcov(nn,ilev) + (I_Zdata(ii,ilev)*this%basis(ii,nn)*this%area(ii))
            end do
         end do
      end do

      ! Multiply by map to get the amplitudes
      !-------------------------------------------
      do ilev=1,nlev
         do nn=1,this%nbas
            O_Bamp(nn,ilev) = 0._r8
            do n2=1,this%nbas
               O_Bamp(nn,ilev) = O_Bamp(nn,ilev) + this%map(n2,nn)*Gcov(n2,ilev)
            end do
         end do
      end do
      deallocate(Gcov)

    end subroutine calc_ZonalProfile_2Damps
    !=======================================================================


    !=======================================================================
    subroutine eval_ZonalProfile_1Dgrid(this,I_Bamp,O_Zdata)
      !
      ! eval_ZonalProfile_1Dgrid: Given the zonal profile basis amplitudes,
      !                           compute 1D data values for the nlat gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalProfile_t):: this
      real(r8),intent(in ):: I_Bamp (:)
      real(r8),intent(out):: O_Zdata(:)
      !
      ! Local Values
      !--------------
      integer:: ii,nn

      ! Construct grid values from basis amplitudes.
      !--------------------------------------------------
      O_Zdata(1:this%nlat) = 0._r8
      do nn=1,this%nbas
        do ii=1,this%nlat
          O_Zdata(ii) = O_Zdata(ii) + (I_Bamp(nn)*this%basis(ii,nn))
        end do
      end do

    end subroutine eval_ZonalProfile_1Dgrid
    !=======================================================================


    !=======================================================================
    subroutine eval_ZonalProfile_2Dgrid(this,I_Bamp,O_Zdata)
      !
      ! eval_ZonalProfile_2Dgrid: Given the zonal profile basis amplitudes,
      !                           compute 2D data values for the nlat,nlev gridpoints.
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalProfile_t):: this
      real(r8),intent(in ):: I_Bamp (:,:)
      real(r8),intent(out):: O_Zdata(:,:)
      !
      ! Local Values
      !--------------
      integer:: ii,nn,ilev

      integer :: nlev

      nlev = size(I_Bamp,dim=2)

      ! Construct grid values from basis amplitudes.
      !--------------------------------------------------
      O_Zdata(1:this%nlat,1:nlev) = 0._r8
      do nn=1,this%nbas
        do ilev=1,nlev
          do ii=1,this%nlat
            O_Zdata(ii,ilev) = O_Zdata(ii,ilev) + (I_Bamp(nn,ilev)*this%basis(ii,nn))
          end do
        end do
      end do

    end subroutine eval_ZonalProfile_2Dgrid
    !=======================================================================


    !=======================================================================
    subroutine init_ZonalAverage(this,IO_lats,IO_area,I_nlat,GEN_GAUSSLATS,USE_LINEARWGTS)
      !
      ! init_ZonalAverage: Initialize the ZonalAverage data structure for the
      !                    given nlat gridpoints. It is assumed that the domain
      !                    of these gridpoints of the profile span latitudes
      !                    from SP to NP.
      !
      !                    The optional GEN_GAUSSLATS flag allows for the
      !                    generation of Gaussian latitudes. The generated grid
      !                    over-writes the values of IO_lats/IO_area passed by
      !                    the user.
      !
      !                    The optional USE_LINEARWGTS flag allows for the 
      !                    linear weighting of gridpoint values as they are 
      !                    distributed to latitude bins, or to assign grid point 
      !                    values to a single bin with equal weighting.
      !                      (default=TRUE for linear weighting)
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalAverage_t) :: this
      real(r8)     ,intent(inout):: IO_lats(:)
      real(r8)     ,intent(inout):: IO_area(:)
      integer         ,intent(in):: I_nlat
      logical,optional,intent(in):: GEN_GAUSSLATS
      logical,optional,intent(in):: USE_LINEARWGTS
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Clats (:)
      real(r8),allocatable:: Glats (:,:)
      real(r8),allocatable:: BinLat(:)
      real(r8),allocatable:: Asum  (:,:)
      real(r8),allocatable:: Anorm (:)
      real(r8):: area(pcols),rlat
      integer :: nn,jj,ierr, astat
      integer :: ncols,lchnk,cc,jlat
      integer :: nlcols, count
      logical :: generate_lats
      character(len=*), parameter :: subname = 'init_ZonalAverage'

      generate_lats = .false.

      if(present(GEN_GAUSSLATS)) then
        generate_lats = GEN_GAUSSLATS
      endif

      this%LINEARWGTS = .true.
      if(present(USE_LINEARWGTS)) then
        this%LINEARWGTS = USE_LINEARWGTS
      endif

      nlcols = get_nlcols_p()

      ! Allocate space
      !-----------------
      if(allocated(this%area   )) deallocate(this%area)
      if(allocated(this%a_norm )) deallocate(this%a_norm)
      if(allocated(this%area_g )) deallocate(this%area_g)
      if(allocated(this%idx_map)) deallocate(this%idx_map)
      if(allocated(this%idx_wgt)) deallocate(this%idx_wgt)

      this%nlat = I_nlat
      allocate(this%area   (I_nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'this%area')
      allocate(this%a_norm (I_nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'this%a_norm')
      allocate(this%area_g (pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'this%area_g')
      allocate(this%idx_map(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'this%idx_map')
      allocate(this%idx_wgt(pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'this%idx_wgt')

      allocate(Clats (I_nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'Clats')
      allocate(BinLat(I_nlat+1), stat=astat)
      call handle_allocate_error(astat, subname, 'BinLat')
      allocate(Glats (pcols,begchunk:endchunk), stat=astat)
      call handle_allocate_error(astat, subname, 'Glats')
      allocate(Asum  (nlcols,I_nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'Asum')
      allocate(Anorm (I_nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'Anorm')

      ! Optionally create the Latitude Gridpoints
      ! and their associated area weights. Otherwise
      ! they need to be supplied by the user.
      !-----------------------------------------------
      if(generate_lats) then

        ! Create a Gaussin grid from SP to NP
        !--------------------------------------
        call sh_create_gaus_grid(this%nlat,Clats,IO_area,ierr)
        if (ierr/=0) then
           call endrun('init_ZonalAverage: Error creating Gaussian grid')
        end if

        ! Convert generated colatitudes SP->NP to Lats and convert
        ! to degrees and scale the area for global 2D integrals
        !-----------------------------------------------------------
        do nn=1,this%nlat
          IO_lats(nn) = (45._r8*Clats(nn)/qrtrPI) - 90._r8
          IO_area(nn) = IO_area(nn)*twoPI
        end do
      else
        ! Convert Latitudes to SP->NP colatitudes in radians
        !----------------------------------------------------
        do nn=1,this%nlat
          Clats(nn) = (IO_lats(nn) + 90._r8)*qrtrPI/45._r8
        end do
      endif

      ! Copy the Lat grid area weights to the data structure
      !-----------------------------------------------------
      this%area(1:this%nlat) = IO_area(1:this%nlat)

      ! Save a copy of the area weights for each 2D gridpoint
      ! and convert Latitudes to SP->NP colatitudes in radians
      !-------------------------------------------------------
      do lchnk=begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        call get_wght_all_p(lchnk, ncols, area)
        do cc = 1,ncols
          rlat=get_rlat_p(lchnk,cc)
          this%area_g(cc,lchnk) = area(cc)
          Glats      (cc,lchnk) = rlat + halfPI
        end do
      end do

      ! Set boundaries for Latitude bins
      !-----------------------------------
      BinLat(1)           = 0._r8
      BinLat(this%nlat+1) = pi
      do nn=2,this%nlat
        BinLat(nn) = (Clats(nn-1)+Clats(nn))/2._r8
      end do

      if(this%LINEARWGTS) then
        ! Each Gridpoint value is distributed with LINEAR weighting 
        ! to the two neighboring Latitude bins
        !------------------------------------------------------
        do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do cc = 1,ncols
            if(Glats(cc,lchnk).lt.Clats(1)) then
              this%idx_map(cc,lchnk) = 1
              this%idx_wgt(cc,lchnk) = (Glats(cc,lchnk)-Clats(2))/(Clats(2)-Clats(1))
            elseif(Glats(cc,lchnk).gt.Clats(I_nlat)) then
              this%idx_map(cc,lchnk) = I_nlat-1
              this%idx_wgt(cc,lchnk) = (Glats(cc,lchnk)-Clats(I_nlat-1))/(Clats(I_nlat)-Clats(I_nlat-1))
            else
              do jj=1,(I_nlat-1)
                if((Glats(cc,lchnk).ge.Clats(jj  )).and. &
                   (Glats(cc,lchnk).le.Clats(jj+1))      ) then
                  this%idx_map(cc,lchnk) = jj
                  this%idx_wgt(cc,lchnk) = (Clats(jj+1)-Glats(cc,lchnk))/(Clats(jj+1)-Clats(jj))
                  exit
                endif
              end do
            endif
          end do
        end do
      else
        ! Each Gridpoint value is assigned to its Latitude bin
        ! with equal weighting.
        !------------------------------------------------------
        this%idx_wgt(:,:) = 1.0_r8

        ! Loop over 2D gridpoints and determine its lat bin index
        !---------------------------------------------------------
        do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do cc = 1,ncols
            jlat = -1
            if((Glats(cc,lchnk)<=BinLat(2)).and. &
               (Glats(cc,lchnk)>=BinLat(1))      ) then
              jlat = 1
            elseif((Glats(cc,lchnk)>=BinLat(this%nlat)  ).and. &
                   (Glats(cc,lchnk)<=BinLat(this%nlat+1))      ) then
              jlat = this%nlat
            else
              do jj=2,(this%nlat-1)
                if((Glats(cc,lchnk)>BinLat(jj  )).and. &
                   (Glats(cc,lchnk)<=BinLat(jj+1))      ) then
                  jlat = jj
                  exit
                endif
              end do
            endif
            if (jlat<1) then
              call endrun('ZonalAverage init ERROR: jlat not in range')
            endif
            this%idx_map(cc,lchnk) = jlat
          end do
        end do
  
        ! Initialize 2D Area sums for each bin
        !--------------------------------------
        Asum(:,:) = 0._r8
        Anorm(:) = 0._r8
        count = 0
        do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do cc = 1,ncols
            jlat = this%idx_map(cc,lchnk)
            count=count+1
            Asum(count,jlat) = this%area_g(cc,lchnk)
          end do
        end do
  
        call shr_reprosum_calc(Asum, Anorm, count, nlcols, I_nlat, gbl_count=ngcols_p, commid=mpicom)
  
        this%a_norm = Anorm
  
        if(.not.all(Anorm(:)>0._r8)) then
           write(iulog,*) 'init_ZonalAverage -- ERROR in Anorm values: '
           do jlat = 1,I_nlat
              if (.not.Anorm(jlat)>0._r8) then
                 write(iulog,*) ' Anorm(',jlat,'): ', Anorm(jlat)
              endif
           end do
           call endrun('init_ZonalAverage -- ERROR in Anorm values')
        endif
      endif

      ! End Routine
      !------------
      deallocate(Clats)
      deallocate(BinLat)
      deallocate(Glats)
      deallocate(Asum)
      deallocate(Anorm)

    end subroutine init_ZonalAverage
    !=======================================================================


    !=======================================================================
    subroutine final_ZonalAverage(this)
      class(ZonalAverage_t) :: this

      if(allocated(this%area   )) deallocate(this%area)
      if(allocated(this%a_norm )) deallocate(this%a_norm)
      if(allocated(this%area_g )) deallocate(this%area_g)
      if(allocated(this%idx_map)) deallocate(this%idx_map)

    end subroutine final_ZonalAverage
    !=======================================================================


    !=======================================================================
    subroutine calc_ZonalAverage_2DbinAvg(this,I_Gdata,O_Zdata)
      !
      ! calc_ZonalAverage_2DbinAvg: Given 2D data values for ncol gridpoints,
      !                             compute the nlat area weighted binAvg profile
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalAverage_t):: this
      real(r8),intent(in ):: I_Gdata(pcols,begchunk:endchunk)
      real(r8),intent(out):: O_Zdata(:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Asum (:,:)
      integer:: nn,ncols,lchnk,cc,jlat
      integer :: nlcols, count, astat
      character(len=*), parameter :: subname = 'calc_ZonalAverage_2DbinAvg'

      nlcols = get_nlcols_p()


      ! Initialize Zonal profile
      !---------------------------
      allocate(Asum(nlcols,this%nlat), stat=astat)
      call handle_allocate_error(astat, subname, 'Asum')
      Asum(:,:) = 0._r8

      O_Zdata(1:this%nlat) = 0._r8

      ! Compute area-weighted sums
      !-----------------------------
      count = 0
      do lchnk=begchunk,endchunk
        ncols = get_ncols_p(lchnk)
        do cc = 1,ncols
          jlat = this%idx_map(cc,lchnk)
          count=count+1
          Asum(count,jlat) = I_Gdata(cc,lchnk)*this%area_g(cc,lchnk)
        end do
      end do

      call shr_reprosum_calc(Asum,O_Zdata,count, nlcols, this%nlat,gbl_count=ngcols_p, commid=mpicom)

      ! Divide by area norm to get the averages
      !-----------------------------------------
      do nn=1,this%nlat
        O_Zdata(nn) = O_Zdata(nn)/this%a_norm(nn)
      end do

      deallocate(Asum)

    end subroutine calc_ZonalAverage_2DbinAvg
    !=======================================================================


    !=======================================================================
    subroutine calc_ZonalAverage_3DbinAvg(this,I_Gdata,O_Zdata)
      !
      ! calc_ZonalAverage_3DbinAvg: Given 3D data values for ncol,nlev gridpoints,
      !                             compute the nlat,nlev area weighted binAvg profile
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalAverage_t):: this
      real(r8),intent(in ):: I_Gdata(:,:,begchunk:)
      real(r8),intent(out):: O_Zdata(:,:)
      !
      ! Local Values
      !--------------
      real(r8),allocatable:: Gsum(:)
      real(r8),allocatable:: Asum(:,:)
      integer:: nn,ncols,lchnk,cc,jlat
      integer:: Nsum,ilev,ns

      integer :: nlev
      integer :: nlcols, count, astat
      character(len=*), parameter :: subname = 'calc_ZonalAverage_3DbinAvg'

      nlev = size(I_Gdata,dim=2)
      nlcols = get_nlcols_p()

      ! Initialize Zonal profile
      !---------------------------
      Nsum = this%nlat*nlev
      allocate(Gsum(Nsum), stat=astat)
      call handle_allocate_error(astat, subname, 'Gsum')
      allocate(Asum(nlcols,Nsum), stat=astat)
      call handle_allocate_error(astat, subname, 'Asum')
      Asum(:,:) = 0._r8

      O_Zdata(1:this%nlat,1:nlev) = 0._r8

      ! Compute area-weighted sums
      !-----------------------------
      do ilev = 1,nlev
         count = 0
         do lchnk=begchunk,endchunk
            ncols = get_ncols_p(lchnk)
            do cc = 1,ncols
               jlat = this%idx_map(cc,lchnk)
               ns = jlat + (ilev-1)*this%nlat
               count=count+1
               Asum(count,ns) = I_Gdata(cc,ilev,lchnk)*this%area_g(cc,lchnk)
            end do
         end do
      end do

      call shr_reprosum_calc(Asum,Gsum, count, nlcols, Nsum, gbl_count=ngcols_p, commid=mpicom)

      ! Divide by area norm to get the averages
      !-----------------------------------------
      do ilev = 1,nlev
         do nn = 1,this%nlat
            ns = nn + (ilev-1)*this%nlat
            O_Zdata(nn,ilev) = Gsum(ns)/this%a_norm(nn)
         end do
      end do

      deallocate(Gsum)
      deallocate(Asum)

    end subroutine calc_ZonalAverage_3DbinAvg
    !=======================================================================


    !=======================================================================
    subroutine set_ZonalAverage_2Dgrid(this,I_Zdata,O_Gdata)
      !
      ! set_ZonalAverage_2Dgrid: Set the 2D data values given 
      !                          the 1D zonal Average profile values,
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalAverage_t) :: this
      real(r8),intent(in ):: I_Zdata(:)
      real(r8),intent(out):: O_Gdata(pcols,begchunk:endchunk)
      !
      ! Local Values
      !--------------
      integer:: ncols,lchnk,cc

      O_Gdata(:,:) = 0._r8

      ! Construct 2D values from bin averages.
      !----------------------------------------
      if(this%LINEARWGTS) then
        ! Use linear weights for neighboring bins
        !----------------------------------------
        do lchnk=begchunk,endchunk
           ncols = get_ncols_p(lchnk)
           do cc = 1,ncols
             if(this%idx_wgt(cc,lchnk).lt.-1._r8) then
               ! South Pole point
               !------------------
               O_Gdata(cc,lchnk) =       -this%idx_wgt(cc,lchnk) *I_Zdata(this%idx_map(cc,lchnk)  ) &
                                  +(1._r8+this%idx_wgt(cc,lchnk))*I_Zdata(this%idx_map(cc,lchnk)+1)
             elseif(this%idx_wgt(cc,lchnk).gt.1._r8) then
               ! North Pole point
               !------------------
               O_Gdata(cc,lchnk) =        this%idx_wgt(cc,lchnk) *I_Zdata(this%idx_map(cc,lchnk)+1) &
                                  +(1._r8-this%idx_wgt(cc,lchnk))*I_Zdata(this%idx_map(cc,lchnk)  )
             else
               O_Gdata(cc,lchnk) =        this%idx_wgt(cc,lchnk) *I_Zdata(this%idx_map(cc,lchnk)  ) &
                                  +(1._r8-this%idx_wgt(cc,lchnk))*I_Zdata(this%idx_map(cc,lchnk)+1)
             endif
           end do
        end do
      else
        ! Set 2D gridpoint values in each bin equal to the given value
        !---------------------------------------------------------------
        do lchnk=begchunk,endchunk
           ncols = get_ncols_p(lchnk)
           do cc = 1,ncols
             O_Gdata(cc,lchnk) = I_Zdata(this%idx_map(cc,lchnk))
           end do
        end do
      endif

    end subroutine set_ZonalAverage_2Dgrid
    !=======================================================================


    !=======================================================================
    subroutine set_ZonalAverage_3Dgrid(this,I_Zdata,O_Gdata)
      !
      ! set_ZonalAverage_3Dgrid: Set the 3D data values given 
      !                          the 2D zonal Average profile values,
      !=====================================================================
      !
      ! Passed Variables
      !------------------
      class(ZonalAverage_t) :: this
      real(r8),intent(in ):: I_Zdata(:,:)
      real(r8),intent(out):: O_Gdata(:,:,begchunk:)
      !
      ! Local Values
      !--------------
      integer:: ncols,lchnk,cc
      integer:: ll

      integer :: nlev
      nlev = size(O_Gdata,dim=2)

      O_Gdata(:,:,:) = 0._r8

      ! Construct 3D values from bin averages.
      !----------------------------------------
      if(this%LINEARWGTS) then
        ! Use linear weights for neighboring bins
        !----------------------------------------
        do ll = 1,nlev
          do lchnk=begchunk,endchunk
             ncols = get_ncols_p(lchnk)
             do cc = 1,ncols
               if(this%idx_wgt(cc,lchnk).lt.-1._r8) then
                 ! South Pole point
                 !------------------
                 O_Gdata(cc,ll,lchnk) =       -this%idx_wgt(cc,lchnk) *I_Zdata(this%idx_map(cc,lchnk)  ,ll) &
                                       +(1._r8+this%idx_wgt(cc,lchnk))*I_Zdata(this%idx_map(cc,lchnk)+1,ll)
               elseif(this%idx_wgt(cc,lchnk).gt.1._r8) then
                 ! North Pole point
                 !------------------
                 O_Gdata(cc,ll,lchnk) =        this%idx_wgt(cc,lchnk) *I_Zdata(this%idx_map(cc,lchnk)+1,ll) &
                                       +(1._r8-this%idx_wgt(cc,lchnk))*I_Zdata(this%idx_map(cc,lchnk)  ,ll)
               else
                 O_Gdata(cc,ll,lchnk) =        this%idx_wgt(cc,lchnk) *I_Zdata(this%idx_map(cc,lchnk)  ,ll) &
                                       +(1._r8-this%idx_wgt(cc,lchnk))*I_Zdata(this%idx_map(cc,lchnk)+1,ll)
               endif
             end do
          end do
        end do
      else
        ! Set 2D gridpoint values in each bin equal to the given value
        !---------------------------------------------------------------
        do ll = 1,nlev
          do lchnk=begchunk,endchunk
             ncols = get_ncols_p(lchnk)
             do cc = 1,ncols
               O_Gdata(cc,ll,lchnk) = I_Zdata(this%idx_map(cc,lchnk),ll)
             end do
          end do
        end do
      endif

    end subroutine set_ZonalAverage_3Dgrid
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

end module ug_spectralmethods_mod
