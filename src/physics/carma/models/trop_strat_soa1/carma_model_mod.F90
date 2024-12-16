!! This CARMA model is for dust aerosols and is based upon Su & Toon, JGR, 2009;
!! Su & Toon, ACP 2011.
!!
!! These dust are not currently radiatively active and do not replace the dust
!! in CAM; however, this is something that could be done in the future.
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!
!! and adds some local functions used to do sea salt emission:
!!
!!   - CARMA_SurfaceWind()
!!   - WeibullWind()
!!
!! @version April-2020
!! @author  Simone Tilmes, Lin Su, Pengfei Yu, Chuck Bardeen
!!  changes to pervious version: rename PURSULF to PRSULF to be easier read in in CAM
!!  Simone Tilmes Aug5 2023: add Ilaria's diagnostic changes

module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use carma_flags_mod
  use carma_model_flags_mod

  use spmd_utils,     only: masterproc
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_get_field, pbuf_get_index
  use time_manager,   only: is_first_step
  use cam_logfile,    only: iulog

  implicit none

  private

  ! Declare the public methods.
  public CARMAMODEL_CalculateCloudborneDiagnostics
  public CARMAMODEL_CreateOpticsFile
  public CARMAMODEL_DefineModel
  public CARMAMODEL_Detrain
  public CARMAMODEL_DiagnoseBins
  public CARMAMODEL_DiagnoseBulk
  public CARMAMODEL_EmitParticle
  public CARMAMODEL_InitializeModel
  public CARMAMODEL_InitializeParticle
  public CARMAMODEL_OutputBudgetDiagnostics
  public CARMAMODEL_OutputCloudborneDiagnostics
  public CARMAMODEL_OutputDiagnostics
  public CARMAMODEL_WetDeposition

  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 2               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 7               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 20              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 2               !! Number of gases

  ! NOTE: This is for now, when Pengfei has only defined sulfates at one weight percent. In the future,
  ! we may want to expand this to match NMIE_WTP and/or NMIE_RH
  integer, public, parameter      :: NREFIDX  = 1               !! Number of refractive indices per element

  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter :: NMIE_RH  = 10              !! Number of relative humidities for mie calculations
  real(kind=f), public, parameter :: mie_rh(NMIE_RH) = (/ 0.1_f, 0.3_f, 0.5_f, 0.7_f, 0.8_f, 0.85_f, &
                                                          0.9_f, 0.92_f, 0.93_f, 0.95_f /)
  integer, public, parameter :: NMIE_WTP = 13              !! Number of weight percents for mie calculations
  real(kind=f), public , parameter :: mie_wtp(NMIE_WTP) = (/ 0.1_f, 0.3_f, 0.5_f, 0.7_f, 0.8_f, 0.83_f, &
                                                            0.86_f, 0.9_f, 0.92_f, 0.94_f, 0.96_f, 0.98_f, 1._f/)

  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_H2SO4          = 1       !! H2SO4 coposition
  integer, public, parameter      :: I_OC             = 2       !! OC composition
  integer, public, parameter      :: I_SOA            = 3       !! SOA composition
  integer, public, parameter      :: I_BC             = 4       !! BC composition
  integer, public, parameter      :: I_DUST           = 5       !! dust composition
  integer, public, parameter      :: I_SALT           = 6       !! sea salt composition

  integer, public, parameter      :: I_GRP_PRSUL     = 1        !! sulfate aerosol
  integer, public, parameter      :: I_GRP_MXAER     = 2        !! mixed aerosol

  integer, public, parameter      :: I_ELEM_PRSUL     = 1       !! sulfate aerosol;  nameing needs to only have 2 charaters  before the element name to work with
                                                                !! partsof the code reading different elements
  integer, public, parameter      :: I_ELEM_MXAER     = 2       !! aerosol
  integer, public, parameter      :: I_ELEM_MXOC      = 3       !! organics aerosol
  integer, public, parameter      :: I_ELEM_MXSOA     = 4       !! secondary organic aerosol
  integer, public, parameter      :: I_ELEM_MXBC      = 5       !! black carbon
  integer, public, parameter      :: I_ELEM_MXDUST    = 6       !! dust aerosol
  integer, public, parameter      :: I_ELEM_MXSALT    = 7       !! sea salt aerosol

  integer, public, parameter      :: I_GAS_H2O        = 1       !! water vapor
  integer, public, parameter      :: I_GAS_H2SO4      = 2       !! sulphuric acid

  real(kind=f), public, parameter         :: Kappa_OC = 0.5_f      !! hygroscopicity of OC
  real(kind=f), public, parameter         :: Kappa_SOA = 0.5_f     !! hygroscopicity of SOA
  real(kind=f), public, parameter         :: Kappa_BC = 0.1_f
  real(kind=f), public, parameter         :: Kappa_DUST = 0.2_f
  real(kind=f), public, parameter         :: Kappa_SALT = 1.0_f
  real(kind=f), public, parameter         :: Kappa_SULF = 0.5_f

  real(kind=f), public, parameter         :: RHO_obc  = 1.35_f          !! dry density of smoke aerosol
  real(kind=f), public, parameter         :: RHO_DUST = 2.65_f          !! dry density of dust particles (g/cm^3) -Lin Su
  real(kind=f), public, parameter         :: RHO_SALT = 2.65_f          !! dry density of sea salt particles (g/cm)
  real(kind=f), public, parameter         :: RHO_SULFATE  = 1.923_f     !! dry density of sulfate particles (g/cm3)

 ! see CARMA_SmokeEmissionRead
! real(kind=f), allocatable, dimension(:,:)     ::   Chla                                       ! Chlorophy11 data (mg/m3)
  real(r8), allocatable, dimension(:,:,:)       ::   BCnew                              ! #/cm2/s
  real(r8), allocatable, dimension(:,:,:)       ::   OCnew


  ! for sea salt flux calculation
  real(r8), parameter             :: uth_salt = 4._r8                !! threshold wind velocity


  ! for dust calculation
  real(kind=f), parameter         :: rClay = 1e-4_f         !! silt/clay particle radius boundary (cm)

  integer                         :: nClay                  !! Number of clay bins (r < 1 um)
  integer                         :: nSilt                  !! Number of silt bins
  real(kind=f)                    :: clay_mf(NBIN)=-huge(1._f) !! clay mass fraction (fraction)
  real(kind=f), allocatable, dimension(:,:) :: soil_factor  !! Soil Erosion Factor (fraction)
  real(kind=f), public, parameter :: WTMOL_H2SO4    = 98.078479_f    !! molecular weight of sulphuric acid

! NOTE: The WeibullK distribution is not currently supported, since the coefficients are not
! generated. This can be added later.
!  real(r8), allocatable, dimension(:,:) :: Weibull_k            ! Weibull K(nlat,nlon
  real(kind=f), public, parameter     :: rmin_PRSUL     = 3.43e-8_f  ! minimum radius (cm)
  real(kind=f), public, parameter     :: vmrat_PRSUL    = 3.67_f     ! volume ratio
  real(kind=f), public, parameter     :: rmin_MXAER     = 5e-6_f     ! minimum radius (cm)
  real(kind=f), public, parameter     :: vmrat_MXAER    = 2.2588_f    !2.4610_f        ! volume ratio

! Physics buffer index for sulfate surface area density
  integer      :: ipbuf4soa(NBIN) = -1
  integer      :: ipbuf4soacm(NBIN) = -1
  integer      :: ipbuf4soapt(NBIN) = -1
  integer      :: ipbuf4jno2 = -1
  real(kind=f) :: aeronet_fraction(NBIN)  !! fraction of BC dV/dlnr in each bin (100%)
  real(kind=f) :: so4inj_dist(NBIN)       !! SO4 injection distribution across bins using a log normal distr. using r=0.95 and sigma =1.5
  real(kind=f) :: so4inj_dist1(NBIN)      !! SO4 injection distribution across bins using a log normal distr. using r=0.95 and sigma =1.5

  integer :: bc_srfemis_ndx=-1, oc_srfemis_ndx=-1
  integer :: so4_elevemis_ndx=-1
  integer :: carma_dustmap(NBIN)        !! mapping of the CARMA dust bins to the surface dust bins.

  ! define refractive indices dependon composition and wavelength
  !
  ! NOTE: It would be better to read this out of files, but this is how Pengfei set it up, so we
  ! will use this for now.
  !
  ! NOTE: Rather than using the values from Pengfei for the sulfate, use the values from MAM. They
  ! have more precision and differ in the imaginary part below 2 um where Pengfei's are truncated at 0.
  ! The MAM values are consistent with OPAC and truncate at 1e-8.
  !real(kind=f), public :: shellreal(NWAVE)    = (/1.890_f,1.913_f,1.932_f,1.568_f,1.678_f,1.758_f,1.855_f,1.597_f,1.147_f,1.261_f,&
  !                1.424_f,1.352_f,1.379_f,1.385_f,1.385_f,1.367_f,&
  !            1.367_f,1.315_f,1.358_f,1.380_f,1.393_f,1.405_f,1.412_f,1.422_f,1.428_f,1.430_f,&
  !            1.422_f,1.468_f,1.484_f,1.164_f/)
  !
  !real(kind=f), public :: shellimag(NWAVE)    = (/0.220_f,0.152_f,0.085_f,0.223_f,0.195_f,0.441_f,0.696_f,0.695_f,0.459_f,0.161_f,&
  !                0.172_f,0.144_f,0.120_f,0.122_f,0.126_f,0.158_f,&
  !            0.158_f,0.057_f,0.003_f,0.001_f,0.001_f,0.000_f,0.000_f,0.000_f,0.000_f,0.000_f,&
  !            0.000_f,0.000_f,0.000_f,0.551_f/)

  real(kind=f), public, parameter :: shellreal(NWAVE)    = (/ 1.89_f, 1.912857_f, 1.932063_f, 1.586032_f, &
               1.677979_f, 1.757825_f, 1.855336_f, 1.596767_f, 1.146559_f, 1.261314_f, 1.424219_f, &
               1.351645_f, 1.378697_f, 1.385_f, 1.385_f, 1.366909_f, 1.366909_f, 1.314577_f, &
               1.357978_f, 1.380309_f, 1.392645_f, 1.404506_f, 1.412181_f, 1.421632_f, &
               1.427968_f, 1.430335_f, 1.441641_f, 1.467642_f, 1.484_f, 1.164128_f /)

  real(kind=f), public, parameter :: shellimag(NWAVE)    = (/ 0.22_f, 0.15185711_f, 0.08457167_f, 0.22250789_f, 0.19499999_f, &
              0.44068847_f, 0.69594361_f, 0.69466153_f, 0.45876573_f, 0.16060575_f, &
              0.1715766_f , 0.14352135_f, 0.12025213_f, 0.12222873_f, 0.12581848_f, 0.15793008_f, &
              1.57930076e-01_f, 5.66869128e-02_f, 2.88634387e-03_f, 1.49071286e-03_f, &
              5.30385233e-04_f, 1.02977119e-04_f, 1.61967358e-05_f, 1.75122678e-06_f, &
              2.21435655e-08_f, 9.99999994e-09_f, 9.99999994e-09_f, 9.99999994e-09_f, &
              9.99999994e-09_f, 5.51133746e-01_f /)

  real(kind=f), public, parameter :: corerealdst(NWAVE)  = &
             (/2.340_f,2.904_f,1.748_f,1.508_f,1.911_f,1.822_f,2.917_f,1.557_f,1.242_f,1.447_f,&
               1.432_f,1.473_f,1.495_f,1.500_f,1.500_f,1.510_f,&
               1.510_f,1.520_f,1.523_f,1.529_f,1.530_f,1.530_f,1.530_f,1.530_f,1.530_f,1.530_f,&
               1.530_f,1.530_f,1.530_f,1.180_f/)

  real(kind=f), public, parameter :: corerealbc (NWAVE)  = &
            (/2.690_f,2.501_f,2.398_f,2.332_f,2.287_f,2.234_f,2.198_f,2.166_f,2.114_f,2.054_f,&
              2.028_f,1.977_f,1.948_f,1.933_f,1.921_f,1.877_f,&
              1.877_f,1.832_f,1.813_f,1.802_f,1.791_f,1.768_f,1.761_f,1.760_f,1.750_f,1.750_f,&
              1.750_f,1.741_f,1.620_f,2.124_f/)

  real(kind=f), public, parameter :: coreimagdst(NWAVE)  = &
             (/0.700_f,0.857_f,0.462_f,0.263_f,0.319_f,0.260_f,0.650_f,0.373_f,0.093_f,0.105_f,&
               0.061_f,0.025_f,0.011_f,0.008_f,0.007_f,0.018_f,&
               0.018_f,0.028_f,0.012_f,0.008_f,0.007_f,0.006_f,0.005_f,0.004_f,0.004_f,0.006_f,&
               0.014_f,0.024_f,0.030_f,0.101_f/)

  real(kind=f), public, parameter :: coreimagbc(NWAVE)   = &
            (/1.000_f,0.884_f,0.825_f,0.791_f,0.764_f,0.734_f,0.714_f,0.696_f,0.668_f,0.644_f,&
              0.624_f,0.604_f,0.593_f,0.586_f,0.580_f,0.556_f,&
              0.556_f,0.527_f,0.503_f,0.492_f,0.481_f,0.458_f,0.451_f,0.440_f,0.430_f,0.443_f,&
              0.461_f,0.470_f,0.450_f,0.674_f/)

  real(kind=f), public, parameter :: waterreal(NWAVE)    = &
   (/ 1.532_f, 1.523857_f, 1.420063_f, 1.274308_f, &
      1.161387_f, 1.142222_f, 1.232189_f, 1.266436_f, 1.295687_f, 1.320659_f, 1.341516_f, &
      1.315192_f, 1.330235_f, 1.339058_f, 1.350425_f, 1.408042_f, 1.408042_f, 1.324462_f, &
      1.276726_f, 1.301847_f, 1.312051_f, 1.321301_f, 1.322836_f, 1.326836_f, 1.330968_f, &
      1.33367_f, 1.339547_f, 1.348521_f, 1.362_f, 1.290783_f /)

  real(kind=f), public, parameter :: waterimag(NWAVE)    = &
   (/ 0.336_f, 0.36000001_f, 0.42623809_f, 0.40341724_f, &
      0.32062717_f, 0.11484398_f, 0.04710282_f, 0.03901278_f, 0.03373134_f, 0.03437707_f, &
      0.09216518_f, 0.0121094_f, 0.01314786_f, 0.01013119_f, 0.00486624_f, 0.0142042_f, &
      1.42042044e-02_f, 1.57659209e-01_f, 1.51634401e-03_f, 1.15906247e-03_f, &
      2.35527521e-04_f, 1.71196912e-04_f, 2.43626002e-05_f, 3.12758360e-06_f, &
      3.74323598e-08_f, 1.63841034e-09_f, 2.49434956e-09_f, 1.52413800e-08_f, &
      3.35000010e-08_f, 3.43825518e-02_f /)

  real(r8), parameter :: onethird = 1._r8/3._r8

contains

  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_DefineModel(carma, rc)

    use physics_buffer, only: pbuf_add_field, dtype_r8

    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure


    ! Local variables
    integer                            :: LUNOPRT              ! logical unit number for output
    character(len=2)                   :: outputname,outputbin
    logical                            :: do_print             ! do print output?
    complex(kind=f)                    :: refidx(NWAVE, NREFIDX) ! refractice indices

    integer                            :: igroup,ibin
    character(len=8)                   :: sname                ! short (CAM) name

    ! Default return code.
    rc = RC_OK

    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_DefineModel: CARMA_Get failed.")

      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_soilerosion_file = ', carma_soilerosion_file
      if (do_print) write(LUNOPRT,*) '  carma_seasalt_emis = ', trim(carma_seasalt_emis)
      if (do_print) write(LUNOPRT,*) '  carma_dustemisfactor = ', carma_dustemisfactor
    end if

    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.

    !call CARMAGROUP_Create(carma, I_GRP_PURSUL, "sulfate", rmin_PRSUL, vmrat_PRSUL, I_SPHERE, 1._f, .false., &
    !                       rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
    !                       scavcoef=0.1_f, is_sulfate=.true., shortname="PRSULF", icoreshell=0, &
    !                       refidx = refidx, refidxS = refidx, refidxC = refidx, do_mie=.true.,imiertn=I_MIERTN_TOON1981)
    !if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    call CARMAGROUP_Create(carma, I_GRP_PRSUL, "sulfate", rmin_PRSUL, vmrat_PRSUL, I_SPHERE, 1._f, .false., &
                           rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.false., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, is_sulfate=.true., shortname="PRSUL", do_mie=.true., &
                           imiertn=I_MIERTN_TOON1981, iopticstype = I_OPTICS_SULFATE)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')


    !call CARMAGROUP_Create(carma, I_GRP_MIXAER, "mixed aerosol", rmin_MIXAER, vmrat_MIXAER, I_SPHERE, 1._f, .false., &
    !                       rc, do_wetdep=.true., do_drydep=.true., solfac=0.2_f, &
    !                       scavcoef=0.1_f, shortname="CRMIX", refidx=refidx, &
    !                       refidxS=refidxS, refidxC=refidxC, do_mie=.true., &
    !                       irhswell=I_MIX, irhswcomp=I_SWG_URBAN, icoreshell=1,imiertn=I_MIERTN_TOON1981)
    !if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    call CARMAGROUP_Create(carma, I_GRP_MXAER, "mixed aerosol", rmin_MXAER, vmrat_MXAER, I_SPHERE, 1._f, .false., &
                           rc, do_wetdep=.false., do_drydep=.true., solfac=0.2_f, &
                           scavcoef=0.1_f, shortname="MXAER", irhswell=I_PETTERS, do_mie=.true., imiertn=I_MIERTN_TOON1981, &
                           iopticstype = I_OPTICS_MIXED_YU_H2O, &
                           neutral_volfrc=-1._f)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')


    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    refidx(:,1) = CMPLX(shellreal(:), shellimag(:), kind=f)
    call CARMAELEMENT_Create(carma, I_ELEM_PRSUL, I_GRP_PRSUL, "Sulfate", &
                             RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="PRSULF", refidx=refidx)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXAER,  I_GRP_MXAER, "Sulfate in mixed sulfate", &
                             RHO_SULFATE, I_VOLATILE, I_H2SO4, rc,  kappa=Kappa_SULF, shortname="MXSULF", refidx=refidx)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXOC,   I_GRP_MXAER, "organic carbon", &
                             RHO_obc, I_COREMASS, I_OC, rc, kappa=Kappa_OC, shortname="MXOC")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXSOA,   I_GRP_MXAER, "secondary organic aerosol", &
                             RHO_obc, I_COREMASS, I_SOA, rc, kappa=Kappa_SOA, shortname="MXSOA")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    refidx(:,1) = CMPLX(corerealbc(:), coreimagbc(:), kind=f)
    call CARMAELEMENT_Create(carma, I_ELEM_MXBC,   I_GRP_MXAER, "black carbon", &
                             RHO_obc, I_COREMASS, I_BC, rc, kappa=Kappa_BC, shortname="MXBC", refidx=refidx)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    refidx(:,1) = CMPLX(corerealdst(:), coreimagdst(:), kind=f)
    call CARMAELEMENT_Create(carma, I_ELEM_MXDUST, I_GRP_MXAER, "dust", &
                             RHO_DUST, I_COREMASS, I_DUST, rc,  kappa=Kappa_DUST, shortname="MXDUST", refidx=refidx)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXSALT, I_GRP_MXAER, "SALT in mixed sulfate", &
                             RHO_SALT, I_COREMASS, I_SALT, rc, kappa=Kappa_SALT, shortname="MXSALT")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')


    ! Define the Solutes



    ! Define the Gases
    refidx(:,1) = CMPLX(waterreal(:), waterimag(:), kind=f)
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, &
                         rc, shortname = "Q", ds_threshold=-0.2_f, refidx=refidx)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
                          I_GCOMP_H2SO4, rc, shortname = "H2SO4")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')


    ! Define the Processes

    call CARMA_AddGrowth(carma, I_ELEM_PRSUL, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddGrowth(carma, I_ELEM_MXAER, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_PRSUL, I_ELEM_PRSUL, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_PRSUL, I_GRP_PRSUL, I_COLLEC_FUCHS, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_MXAER, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    !----------------- add pbuf ------------------
    do igroup = 1, NGROUP

      call CARMAGROUP_Get(carma, igroup, rc, shortname=sname)
      if (rc < 0) call endrun('carma_register::CARMAGROUP_Get failed.')
      !write(*,*) "igroup",igroup,"sname",sname

      ! sulfate mass and number density for each bin
      ! e.g. CRSULF01 first element mass mixing ratio; NBMXAER01 #/kg
      do ibin=1,NBIN
         write (outputbin, "(I2.2)") ibin
         if (igroup==I_GRP_MXAER) then
           call pbuf_add_field("DQDT_MXSOA"//outputbin,'global',dtype_r8,(/pcols,pver/), ipbuf4soa(ibin))
           call pbuf_add_field("MXSOA"//outputbin//"CM",'physpkg',dtype_r8,(/pcols,pver/), ipbuf4soacm(ibin))
           call pbuf_add_field("MXSOA"//outputbin//"PT",'physpkg',dtype_r8,(/pcols,pver/), ipbuf4soapt(ibin))
         end if
      end do
   end do

    ! no2 photolysis rate constant (/sec)
    call pbuf_add_field('JNO2', 'global', dtype_r8, (/pcols,pver/), ipbuf4jno2)

    !---------------------------------------------

    return
  end subroutine CARMAMODEL_DefineModel


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009
  !!  @author  Chuck Bardeen
  !!
  !!  @see CARMASTATE_SetDetrain
  subroutine CARMAMODEL_Detrain(carma, cstate, cam_in, dlf, state, icol, dt, rc, rliq, prec_str, snow_str, &
     tnd_qsnow, tnd_nsnow)
    use camsrfexch, only: cam_in_t

    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in           !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(r8), intent(in)                 :: dt               !! time step (s)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s)
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional      :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(out), optional      :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)

    ! Default return code.
    rc = RC_OK

    return
  end subroutine CARMAMODEL_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, rliq, prec_str, snow_str)

    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)       :: state        !! physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s)
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)

    ! local variables
    real(r8), pointer, dimension(:,:)    :: dqdt_soa              !! soa tendency due to gas-aerosol exchange  kg/kg/s
    real(r8), pointer, dimension(:,:)    :: jno2_rate             !! jno2 tendency due to gas-aerosol exchange  kg/kg/s
    real(r8), pointer, dimension(:,:)    :: soacm                 !! aerosol tendency due to gas-aerosol exchange  kg/kg/s
    real(r8), pointer, dimension(:,:)    :: soapt                 !! aerosol tendency due to no2 photolysis  kg/kg/s
    real(r8)                             :: mmr_core(cstate%f_NZ)!! mass mixing ratio of the core (kg/kg)
    real(r8)                             :: mmr_soa(cstate%f_NZ)  !! mass mixing ratio of soa element (kg/kg)
    real(r8)                             :: mmr(cstate%f_NZ)      !! mass mixing ratio per bin (kg/kg)
    real(r8)                             :: delta_soa(cstate%f_NZ)     !! mass mixing ratio differences from soa gas-aerosol-exchange
    integer                              :: icorelem(NELEM), ncore,ienconc,icore, ielem, ielem_soa, igroup, ibin, icomposition, n, err

    ! Default return code.
    rc = RC_OK

    ! get no2 photolysis rates if they exist
    call pbuf_get_field(pbuf, ipbuf4jno2, jno2_rate)     ! surface area density

    ! get SOA tendency pbuf field for the mixed group and every bin

    igroup = I_GRP_MXAER
    call CARMAGROUP_Get(carma, igroup, rc, ienconc=ienconc, ncore=ncore, icorelem=icorelem)
    if (rc /= RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

    do ibin = 1, NBIN

      ! Iterate over the core elements, looking for the SOA element. Once found,
      ! determine the new SOA taking into account both the addition of condensed
      ! SOA and the loss of photolyzed SOA.
      do ielem = 1, ncore

        call CARMASTATE_GetBin(cstate, icorelem(ielem), ibin, mmr(:), rc)
        if (rc /= RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_GetBin failed.')

        call CARMAELEMENT_GET(carma, icorelem(ielem), rc, icomposition=icomposition)
        if (rc /= RC_OK) call endrun('CARMA_DiagnoseBins::CARMAELEMENT_Get failed.')

        ! Only need to make adjustments for the SOA.
        if (icomposition == I_SOA) then
          call pbuf_get_field(pbuf, ipbuf4soa(ibin), dqdt_soa)     ! surface area density

          ! Add that soa tendency from chemistry to the aerosol.
          !
          !   NOTE: dqdt is in kg/kg/s
          mmr_soa(:) = mmr(:) + dqdt_soa(icol,:) * dt

          ! Save the chemistry tendency so it can by output in the diagnostics.
          call pbuf_get_field(pbuf, ipbuf4soacm(ibin), soacm)
          soacm(icol,:) = dqdt_soa(icol,:)

          ! Save the NO2 photolysis tendency so it can by output in the diagnostics.
          !
          ! NOTE: Simone, what is the 0.0004_r8??
          call pbuf_get_field(pbuf, ipbuf4soapt(ibin), soapt)
          soapt(icol,:) = - 0.0004_r8 * jno2_rate(icol,:) * mmr_soa(:)

          ! Now adjust the SOA for the loss by the photolysis rate provided by the
          ! chemistry.
          mmr_soa(:) = max(0.0_r8, mmr_soa(:) + soapt(icol,:) * dt)

          ! Save out these new values for SOA.
          call CARMASTATE_SetBin(cstate, icorelem(ielem), ibin, mmr_soa, rc)
          if (rc /= RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_SetBin failed.')

          exit
        end if  !mxsoa
      end do  !ielem
    end do  !nbin

  end subroutine CARMAMODEL_DiagnoseBins


  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch, only: cam_out_t

    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    type(cam_out_t),      intent(inout)  :: cam_out   !! cam output to surface models
    type(physics_state), intent(in)      :: state     !! physics state variables
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer
    type(physics_ptend), intent(inout)   :: ptend     !! constituent tendencies
    integer, intent(in)                  :: icol      !! column index
    real(r8), intent(in)                 :: dt        !! time step
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s)
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(inout), optional    :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(inout), optional    :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    real(r8), intent(out), optional      :: re_ice(pcols,pver)    !! ice effective radius (m)

   ! Local variables
    real(r8)                             :: numberDensity(cstate%f_NZ)
    real(r8)                             :: totad(cstate%f_NZ)
    real(r8)                             :: ad(cstate%f_NZ)       !! aerosol wet surface area density (cm2/cm3)
    real(r8)                             :: totreff(cstate%f_NZ)  !! total volume density, used to calculate total effective radius (cm) for history output
    real(r8)                             :: reff(cstate%f_NZ)     !! wet effective radius (m)
    real(r8)                             :: mmr(cstate%f_NZ)      !! mass mixing ratio per bin (kg/kg)
    real(r8)                             :: coremmr(cstate%f_NZ)  !! mmr of all the core
    real(r8)                             :: mmr_gas(cstate%f_NZ)  !! gas mass mixing ratio (kg/kg)
    real(r8)                             :: numnkg(cstate%f_NZ)   !! total number density (#/kg)
    real(r8)                             :: r_wet(cstate%f_NZ)    !! Sulfate aerosol bin wet radius (cm)
    real(r8)                             :: elem1mr(cstate%f_NZ)  !! First element mass mixing ratio (kg/kg)
    real(r8)                             :: binnkg(cstate%f_NZ)   !! number density per bin (#/kg)
    real(r8)                             :: kappa(cstate%f_NZ)    !! hygroscopicity parameter (Petters & Kreidenweis, ACP, 2007)
    real(r8)                             :: rhoa_wet(cstate%f_NZ) !! wet air density (kg/m3)
    real(r8)                             :: wtpct(cstate%f_NZ)    !! sulfate weight percent
    real(r8)                             :: rmass(NBIN)           !! dry mass
    real(r8)                             :: rhop_dry(cstate%f_NZ) !! dry particle density [g/cm3]

    integer                              :: ibin, igroup, igas, icomposition
    integer                              :: icorelem(NELEM), ncore,ienconc,icore
    character(len=8)                     :: sname                 !! short (CAM) name

    real(r8), pointer, dimension(:,:)    :: sadsulf_ptr           !! Total surface area density pointer (cm2/cm3)
    real(r8), pointer, dimension(:,:)    :: reffaer_ptr           !! Total effective radius pointer (cm) for history output
    real(r8), pointer, dimension(:,:)    :: wtp_ptr               !! weight percent pointer
    real(r8), pointer, dimension(:,:)    :: sad_ptr               !! Surface area density pointer
    real(r8), pointer, dimension(:,:)    :: reff_ptr              !! Effective radius pointer
    real(r8), pointer, dimension(:,:)    :: numnkg_ptr            !! Each group number density pointer
    real(r8), pointer, dimension(:,:)    :: binnkg_ptr            !! Each bin number density pointer
    real(r8), pointer, dimension(:,:)    :: elem1mr_ptr           !! First element mmr pointer
    real(r8), pointer, dimension(:,:)    :: kappa_ptr             !! kappa pointer
    real(r8), pointer, dimension(:,:)    :: wetr_ptr              !! wet radius pointer
    real(r8), pointer, dimension(:,:)    :: dryr_ptr              !! dry radius

    ! Default return code.
    rc = RC_OK

    return
  end subroutine CARMAMODEL_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Lin Su, Pengfei Yu, Chuck Bardeen
  !! @version Dec-2010
  subroutine CARMAMODEL_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, pbuf, rc)
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_lon_all_p, get_lat_all_p
    use time_manager,  only: get_curr_date, get_perp_date, is_perpetual
    use camsrfexch,    only: cam_in_t
    use cam_history,   only: outfld

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    type(physics_buffer_desc), pointer :: pbuf(:)               !! physics buffer
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer      :: ilat(pcols)             ! latitude index
    integer      :: ilon(pcols)             ! longitude index
    real(r8)     :: clat(pcols)             ! latitude
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: p                       ! plev index
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    real(r8)     :: smoke(pcols)            ! smoke emission flux (molecues/cm2/s)
    real(r8)     :: rhoa(pcols,pver)        ! density of air  g/cm3
    real(r8)     :: so4_inj(pcols,pver)     ! so4 emission flux (molecues/cm3/s)
    real(r8)     :: so4_tendency_factor(pcols,pver)     ! Convertion factor from molec/cm3/s to kg/kg/s
    integer      :: igroup                  ! the index of the carma aerosol group
    character(len=32) :: shortname          ! the shortname of the group



    ! -------- local variables added for dust and sea-salt model ------------
    real(r8)            :: ch                                 ! dimensional factor & tuning number,
    real(r8)            :: rmass(NBIN)                        ! bin mass (g)
    real(r8)            :: r                                  ! bin center (cm)
    real(r8)            :: rdust                              ! dust bin center (cm)
    real(r8)            :: dustFlux                           ! dust flux (kg/m2/s)
    real(r8)            :: rsalt                              ! salt bin center (cm)
    real(r8)            :: drsalt                             ! salt bin width (cm)
    real(r8)            :: rhop(NBIN)                         ! element density (g/cm3)
    real(r8)            :: vrfact
    real(r8)            :: uth                                ! threshold wind velocity (m/s)
    real(r8)            :: uv10                               ! 10 m wind speed (m/s)
    real(r8)            :: cd10                               ! 10-m drag coefficient ()
    real(r8)            :: wwd                                ! raw wind speed (m/s)
    real(r8)            :: sp                                 ! mass fraction for soil factor
    integer             :: idustbin                           ! ibin to use for dust production, smallest silt bin for clay

! ------------ local variables added for organics model ----------------------
    real(r8)     :: dr
    real(r8)     :: aeronet(NBIN)                       ! AERONET DATA, Sep.20, 2002, Jaru Reserve, Brazil (refer to MATICHUK et al., 2008)
    real(r8)     :: saltFlux(pcols)                     ! sea salt flux to calculate marine POA
    integer      :: LUNOPRT                             ! logical unit number for output
    logical      :: do_print                            ! do print output?

    real(r8),parameter :: OMtoOCratio = 1.8_r8           ! Need better names and doc
    real(r8),parameter :: SmoketoSufaceFlux = 1.9934e-22_r8 ! SmoketoSufaceFlux = BC molecular weight
                                                            ! (12 g/mol)/avocadro constant (6e-23 #/mol) *10
    real(r8), pointer :: BCemis_ptr(:), OCemis_ptr(:)
    real(r8), pointer :: SO4elevemis_ptr(:,:)

    ! Default return code.
    rc = RC_OK
    smoke(:) = -huge(1._r8)
    so4_inj(:,:) = -huge(1._r8)
    ch = carma_dustemisfactor

    ! Determine the day of year.
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8

    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.
    tendency(:ncol, :pver) = 0.0_r8

     ! Add Emission (surfaceFlux) here.

    !!*******************************************************************************************************

    !! add an element, first element is total number with emission from both OC and BC;
    !! second element is BC mass
    !! by Pengfei Yu
    !! Feb.22 2012
    !!*******************************************************************************************************


    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
    if (RC < RC_ERROR) return

    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, rmass=rmass)
    if (RC < RC_ERROR) return

     !!*******************************************************************************************************

    !if (masterproc) then
    !  call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    !
    ! if (do_print) then
    !   write(carma%f_LUNOPRT,*) 'AERONET', aeronet
    !   write(carma%f_LUNOPRT,*) 'dr', dr
    !   write(carma%f_LUNOPRT,*) 'r', r
    ! end if
    !end if

    !!*******************************************************************************************************

    if(carma_BCOCemissions == 'Specified')then
      call pbuf_get_field(pbuf, bc_srfemis_ndx, BCemis_ptr)
      call pbuf_get_field(pbuf, oc_srfemis_ndx, OCemis_ptr)
    end if
    if(carma_SO4elevemis== 'Specified')then
      call pbuf_get_field(pbuf, so4_elevemis_ndx, SO4elevemis_ptr)
    end if

    ! Organic carbon emssions
    if (ielem == I_ELEM_MXOC) then
       if (carma_BCOCemissions == 'Yu2015') then
          call get_lat_all_p(lchnk, ncol, ilat)
          call get_lon_all_p(lchnk, ncol, ilon)
          do icol = 1,ncol
             smoke(icol) = OCnew(ilat(icol), ilon(icol), mon)*OMtoOCratio
          end do
       elseif(carma_BCOCemissions == 'Specified')then
          smoke(:ncol) = OCemis_ptr(:ncol)
       end if

!  st  scip Fsub PBAFlux etcfor now
       surfaceFlux(:ncol) = surfaceFlux(:ncol) + smoke(:ncol)*aeronet_fraction(ibin)*SmoketoSufaceFlux
    end if

    ! Black carbon emissions
    if (ielem == I_ELEM_MXBC) then
       if (carma_BCOCemissions == 'Yu2015') then
          do icol = 1,ncol
             smoke(icol) = BCnew(ilat(icol), ilon(icol), mon)
          end do
       elseif(carma_BCOCemissions == 'Specified') then
          smoke(:ncol) = BCemis_ptr(:ncol)
       end if

       surfaceFlux(:ncol) = surfaceFlux(:ncol) + smoke(:ncol)*aeronet_fraction(ibin)*SmoketoSufaceFlux
    end if

    if(carma_SO4elevemis == 'Specified') then
       ! Sulfate emissions
       if (ielem == I_ELEM_PRSUL)  then
          ! convert from #/kg to kg/kg  = 1.e-3 *  mw/avog (6e-23)    !kg/kg
          ! convert from #/cm3/s to kg/kg/s = 1.e3 * density of air * mw / avog
          !AVG: molec/mol R_AIR: units?
          !rhoa
          !number Density
          !rhoa(:ncol,:) = 10._r8 * state%pmid(:ncol,:) / (R_AIR * state%t(:ncol,:))
          !pmid is in Pa (Pa->dynes (factor of 10.), T (K), -> g/cm3

          !so4_tendency_factor(:ncol,:) = rhoa(:ncol,:) * WTMOL_H2SO4 / AVG  !molec/cm3/s to kg/kg

          so4_inj(:ncol,:) = SO4elevemis_ptr(:ncol,:)


          ! set so4_inj larger 0. because of potential negative missing values
          do icol = 1,ncol
             do p = 1,pver
                rhoa(icol,p) = 10._r8 * state%pmid(icol,p) / (R_AIR * state%t(icol,p))
                !pmid is in Pa (Pa->dynes (factor of 10.), T (K), -> g/cm3
                !emis = molec/cm3/s
                !rhoa = g/cm3
                !mw = g/mol
                !avg =  molec/mol
                !so4_tendency_factor(icol,p) = rhoa(icol,p) * WTMOL_H2SO4 / AVG  !molec/cm3/s to kg/kg
                so4_tendency_factor(icol,p) =  WTMOL_H2SO4 / AVG / rhoa(icol,p)  !molec/cm3/s to kg/kg
                so4_inj(icol,p) = max(0._r8,so4_inj(icol,p))
                if (so4_inj(icol,p).gt.0._r8) then
                   tendency(icol,p) = so4_inj(icol,p)*so4inj_dist(ibin)*so4_tendency_factor(icol,p)
                end if
             end do
          end do
       end if
    end if

    ! Dust emissions
    if (ielem == I_ELEM_MXDUST) then

      ! The radius should be determined by the dust density not the group
      ! density
      call CARMAELEMENT_Get(carma, I_ELEM_MXDUST, rc, rho=rhop)
      if (RC < RC_ERROR) return

      ! Calculate the radius assuming that all the mass will be emitted as this
      ! element.
      rdust = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin)) ** (1._r8 / 3._r8)

      ! Is this clay or silt?
      !
      ! NOTE: It is assumed that 90% of the mass will be silt and 10% will
      ! be clay.
      !
      ! NOTE: For clay bins, use the smallest silt bin to calculate the
      ! mass and then scale that into each clay bin based upon interpolation of
      ! Tegen and Lacis [1996].
      if (rdust >= rClay) then
        sp         = 0.9_r8 / nSilt
        idustbin   = ibin
      else
        sp         = 0.1_r8 / nClay
        idustbin   = nClay + 1
      end if

      ! Process each column.
      do icol = 1,ncol

        call CARMAMODEL_SurfaceWind(carma, icol, I_ELEM_MXDUST, igroup, idustbin, cam_in, uv10, wwd, uth, rc)

        ! Is the wind above the threshold for dust production?
        if (sqrt(wwd) > uth) then
          dustFlux = ch * soil_factor(icol, lchnk) * sp * &
                              wwd * (sqrt(wwd) - uth)
        else
          dustFlux = 0._r8
        endif

        ! Scale the clay bins based upon the smallest silt bin.
        dustFlux = clay_mf(ibin) * dustFlux

        ! Add the dust flux to the accumulated emissions (important for I_ELEM_MXAER)
        surfaceFlux(icol) = surfaceFlux(icol) + dustFlux
      end do

      ! For debug purposes, output the soil erosion factor.
      call outfld('CRSLERFC', soil_factor(:ncol, lchnk), ncol, lchnk)
    end if


    ! Sea salt emissions
    if (ielem == I_ELEM_MXSALT) then

      ! The radius should be determined by the dust density not the group
      ! density
      call CARMAELEMENT_Get(carma, I_ELEM_MXSALT, rc, rho=rhop)
      if (RC < RC_ERROR) return

      ! Calculate the radius assuming that all the mass will be emitted as sea
      ! salt.
      vrfact = ((3._r8/2._r8 / PI / (vmrat_MXAER + 1._r8))**(1._r8 / 3._r8)) * ((vmrat_MXAER**(1._r8 / 3._r8)) - 1._r8)
      rsalt = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin))**(1._r8 / 3._r8)
      drsalt = vrfact * ((rmass(ibin)/rhop(ibin))**(1._r8 / 3._r8))

      ! get sea spray aerosol flux first (for ibin; SaltFlux(:ncol) unit:kg/m2/s)
      call CARMAMODEL_SaltFlux(carma, ibin, state, rsalt, drsalt, rmass(ibin), cam_in, saltFlux, rc)

!st  not used currently  but done by Pengfei
       !! introduce marine POA emission, use ChlorophyII-dependent mass contribution of OC
       !! see Gantt et al., 2009
       !! for sub-micron, I use sea salt flux instead of sub-micron marine particles
       !! needed to verify later
       !! Added by Pengfei Yu
       !! Oct.6.2012
       ! get [Chl-a] data
  !!   do icol = 1, ncol
  !!       if (Chla(ilat(icol), ilon(icol)) .lt. 0._r8) then
  !!          Fsub(icol) = 0._r8
  !!       else
  !!          Fsub(icol) = Chla(ilat(icol), ilon(icol)) * 0.63_r8 + 0.1_r8
  !!       endif
  !!       Fsub(icol) = min(Fsub(icol), 1._r8)
  !!   enddo
  !!   surfaceFlux(:ncol) = SaltFlux(:ncol)
  !!   ! sea salt (NaCl) flux should exclude marine organics and marine sulfate
  !!   if (carma%f_group(igroup)%f_r(ibin) .le. 0.5e-4_r8) then
  !!       !surfaceFlux(:ncol) = SaltFlux(:ncol)*(1._r8-0.0983_r8) - SaltFlux(:ncol) * Fsub(:ncol)
  !!        surfaceFlux(:ncol) = (SaltFlux(:ncol) - SaltFlux(:ncol)*Fsub(:ncol))/1.0983_r8
  !!   else
  !!       !surfaceFlux(:ncol) = SaltFlux(:ncol)*(1._r8-0.0983_r8) - SaltFlux(:ncol) * (Fsub(:ncol)*0.03_r8)
  !!        surfaceFlux(:ncol) = (SaltFlux(:ncol) - SaltFlux(:ncol)*Fsub(:ncol)*0.03_r8)/1.0983_r8
  !!   endif
      surfaceFlux(:ncol) = surfaceFlux(:ncol) + saltFlux(:ncol)
    end if

    return
  end subroutine CARMAMODEL_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMAMODEL_InitializeModel(carma, lq_carma, pbuf2d, rc)
    use cam_history,  only: addfld,  horiz_only, add_default
    use constituents, only: pcnst

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! -------- local variables ----------
    integer            :: ibin                                ! CARMA bin index
    real(r8)           :: r(NBIN), dr(NBIN), rdust(NBIN),robc(NBIN),drobc(NBIN),rm(NBIN),rhop(NBIN)       ! bin center (cm)
    integer            :: count_Silt                          ! count number for Silt
    integer            :: igroup                              ! the index of the carma aerosol group
    integer            :: ielem                               ! the index of the carma aerosol element
    character(len=32)  :: shortname                           ! the shortname of the element
    integer            :: LUNOPRT                             ! logical unit number for output
    logical            :: do_print                            ! do print output?

    integer :: i, idata,isizebin,ibin_local
    integer,parameter :: aeronet_dim1 = 22
    integer,parameter :: aeronet_dim2 = 4
    real(r8),dimension(aeronet_dim1,aeronet_dim2) :: sizedist_aeronet
    real(r8),dimension(aeronet_dim1) :: sizedist_avg
    real(r8),dimension(NBIN) :: sizedist_carmabin
    real(r8) :: rmass(NBIN) !! dry mass
    real(r8) :: vrfact
    real(r8) :: rgeo
    real(r8) :: siglog, siglogsq, sq2pi
    character(len=16)    :: binname      !! names bins

    real(r8),parameter :: size_aeronet(aeronet_dim1) = (/0.050000_r8,0.065604_r8,0.086077_r8,0.112939_r8,0.148184_r8, &
         0.194429_r8,0.255105_r8,0.334716_r8,0.439173_r8,0.576227_r8,0.756052_r8,0.991996_r8,1.301571_r8,1.707757_r8, &
         2.240702_r8,2.939966_r8,3.857452_r8,5.061260_r8,6.640745_r8,8.713145_r8,11.432287_r8,15.000000_r8/)*1.e-4_r8 !um to cm

    ! Default return code.
    rc = RC_OK

    ! Determine how many clay and how many silt bins there are, based
    ! upon the bin definitions and rClay.
    !
    ! TBD: This should use the radii rather than being hard coded.
    ! nClay = 8
    ! nSilt = NBIN - nClay
    do ielem = 1, NELEM

       ! To get particle radius, need to derive from rmass and density of dust.
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname, rho=rhop)
       if (RC < RC_ERROR) return

       call CARMAGROUP_GET(carma, igroup, rc, rmass=rmass)
       if (RC < RC_ERROR) return

       if (shortname .eq. "MXDUST") then

          count_Silt = 0
          do ibin = 1, NBIN

             ! Calculate the radius assuming that all the mass will be emitted as this
             ! element.
             rdust(ibin) = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin))**(1._r8 / 3._r8)

             if (rdust(ibin) >= rclay) then
                count_Silt = count_Silt + 1
             else
             end if
          end do
          nSilt = count_Silt
          nClay = NBIN - nSilt
       end if
    end do

    ! Read in the soil factors.
    call CARMAMODEL_ReadSoilErosionFactor(rc)
    if (RC < RC_ERROR) return

    ! To determine Clay Mass Fraction
    do ielem = 1, NELEM
       ! To get particle radius
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
       if (RC < RC_ERROR) return

       if (shortname .eq. "MXDUST") then
          call CARMAMODEL_ClayMassFraction(carma, igroup, rdust, rc)
       end if
    end do

    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")

      if (do_print) then
        write(carma%f_LUNOPRT,*) 'Initializing CARMA dust model ...'
        write(carma%f_LUNOPRT,*) 'nClay = ', nClay, ' nSilt = ', nSilt
        write(carma%f_LUNOPRT,*) 'clay_mf = ', clay_mf
        write(carma%f_LUNOPRT,*) 'soil_factor = ', soil_factor

        write(carma%f_LUNOPRT,*) 'CARMA dust initialization complete'
      end if
    end if

    call addfld('CRSLERFC', horiz_only, 'A', 'fraction', 'CARMA soil erosion factor')

    if (carma_BCOCemissions == 'Yu2015')then
       ! Added by Pengfei Yu to read smoke emission data
       call CARMAMODEL_BCOCread(rc)
    end if
    if(carma_BCOCemissions == 'Specified')then
       bc_srfemis_ndx = pbuf_get_index("BC_srfemis")
       oc_srfemis_ndx = pbuf_get_index("OC_srfemis")
    end if

    ! prescribed sulfate emissions for stratospheric aerosol injections
    if(carma_SO4elevemis == 'Specified')then
       so4_elevemis_ndx = pbuf_get_index("SO4_elevemis")
    end if

    if (is_first_step()) then

       ! Initialize physics buffer fields
       do igroup = 1, NGROUP
          do ibin = 1, NBIN
             if (igroup==I_GRP_MXAER) then
                call pbuf_set_field(pbuf2d, ipbuf4soa(ibin), 0.0_r8 )
             end if
          end do
       end do

       call pbuf_set_field(pbuf2d, ipbuf4jno2, 0.0_r8 )
    endif

    sizedist_aeronet(:aeronet_dim1,1) = (/0.000585_r8,0.006080_r8,0.025113_r8,0.052255_r8,0.079131_r8,0.081938_r8, &
         0.035791_r8,0.010982_r8,0.005904_r8,0.007106_r8,0.011088_r8,0.012340_r8,0.010812_r8,0.010423_r8, &
         0.011892_r8,0.016529_r8,0.023967_r8,0.026854_r8,0.017901_r8,0.007226_r8,0.002161_r8,0.000544_r8/)
    sizedist_aeronet(:aeronet_dim1,2) = (/0.000541_r8,0.006524_r8,0.026103_r8,0.050825_r8,0.077730_r8,0.080545_r8, &
         0.035400_r8,0.011143_r8,0.005753_r8,0.006095_r8,0.008730_r8,0.010794_r8,0.011517_r8,0.012051_r8, &
         0.012362_r8,0.014710_r8,0.019738_r8,0.022156_r8,0.014892_r8,0.005976_r8,0.001891_r8,0.000573_r8/)
    sizedist_aeronet(:aeronet_dim1,3) = (/0.000747_r8,0.009291_r8,0.043556_r8,0.099216_r8,0.142377_r8,0.108606_r8, &
         0.043723_r8,0.016385_r8,0.008318_r8,0.005597_r8,0.004431_r8,0.004131_r8,0.004980_r8,0.007484_r8, &
         0.011795_r8,0.017235_r8,0.022404_r8,0.025216_r8,0.022521_r8,0.013752_r8,0.005051_r8,0.001057_r8/)
    sizedist_aeronet(:aeronet_dim1,4) = (/0.000979_r8,0.007724_r8,0.034451_r8,0.090410_r8,0.135893_r8,0.103115_r8, &
         0.046047_r8,0.018989_r8,0.009149_r8,0.005034_r8,0.003199_r8,0.002680_r8,0.003249_r8,0.005105_r8, &
         0.008370_r8,0.012542_r8,0.016973_r8,0.021107_r8,0.022077_r8,0.015639_r8,0.006001_r8,0.001115_r8/)

    sizedist_avg(:) = 0._r8
    do idata = 1,aeronet_dim2
       sizedist_avg(:) = sizedist_avg(:) + sizedist_aeronet(:,idata)
    end do
    sizedist_avg(:) = sizedist_avg(:)*0.25_r8

    do igroup = 1,NGROUP
      call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, rmass=rmass)


      if (shortname .eq. "MXAER") then

        !interpolate into carma bin
        sizedist_carmabin = 0._r8

        do ibin_local = 1, NBIN
          ! Calculate the radius assuming that all the mass will be emitted as this
          ! element.
          vrfact = ((3._r8/2._r8 / PI / (vmrat_MXAER + 1._r8))**(1._r8 / 3._r8)) * ((vmrat_MXAER**(1._r8 / 3._r8)) - 1._r8)
          robc(ibin_local) = (3._r8 * rmass(ibin_local) / 4._r8 / PI / rho_obc)**(1._r8 / 3._r8)
          drobc(ibin_local) = vrfact * ((rmass(ibin_local)/rho_obc) **(1._r8 / 3._r8))

          if(robc(ibin_local) .lt. size_aeronet(1)) then
            sizedist_carmabin(ibin_local) = sizedist_avg(1)
          end if
          if(robc(ibin_local) .ge. size_aeronet(aeronet_dim1)) then
            sizedist_carmabin(ibin_local) = sizedist_avg(aeronet_dim1)
          end if
          do isizebin= 1,aeronet_dim1-1
            if( robc(ibin_local) .ge. size_aeronet(isizebin) .and.  robc(ibin_local) .lt. size_aeronet(isizebin+1))then
              sizedist_carmabin(ibin_local) = sizedist_avg(isizebin)*(size_aeronet(isizebin+1)-robc(ibin_local))/&
                  (size_aeronet(isizebin+1)-size_aeronet(isizebin))&
                  +sizedist_avg(isizebin+1)*(robc(ibin_local)-size_aeronet(isizebin))&
                  /(size_aeronet(isizebin+1)-size_aeronet(isizebin))
            end if
          end do
        end do

        rm(:) = 0._r8
        do ibin_local = 1, NBIN
          rm(ibin_local) = sizedist_carmabin(ibin_local)*drobc(ibin_local)/robc(ibin_local)*RHO_obc*1.e-15_r8         ! kg
        enddo

        do ibin_local = 1, NBIN
          aeronet_fraction(ibin_local) = rm(ibin_local)/sum(rm(:))
        end do

      end if
    end do

    ! Produce lognormal size distribtuion for sulfate emissions (SO4 geoengienering experiments)

    ! Define specific for SO4 injection, e.g.,mean dry radius: 0.095, sigma = 1.5
    so4inj_dist(:) = 0.0_r8
    so4inj_dist1(:) = 0.0_r8
    rgeo=0.095e-4_f                    ! mean radius for aerosol injections in cm
    siglog=log(1.5_r8)    ! assumed log normal distribtuion around mean radius for aerosol injections
    siglogsq=siglog**2_f
    sq2pi = sqrt(2._r8*pi)
    !aer_Vrat =  vmrat_PRSUL

    call CARMAGROUP_GET(carma, I_GRP_PRSUL, rc, r=r, dr=dr, shortname=shortname, rmass=rmass)

    !interpolate into carma bin

    do ibin_local = 1, NBIN
       ! Size Distribution-Parameter: log-normal distribution applied using Seinfeld and Pandis (2016)
       so4inj_dist1(ibin_local)=dr(ibin_local)/(r(ibin_local)*sq2pi*siglog)*exp(-(((log(r(ibin_local)/rgeo))**2._r8)/(2._r8*siglogsq)))
       so4inj_dist(ibin_local)=dr(ibin_local)/(r(ibin_local)*sq2pi*siglog)*exp(-(((log(r(ibin_local)/rgeo))**2._r8)/(2._r8*siglogsq)))
       so4inj_dist1(ibin_local) = so4inj_dist1(ibin_local) *rmass(ibin_local)
    end do
    so4inj_dist(:) = so4inj_dist(:) / sum(so4inj_dist)
    so4inj_dist1(:) = so4inj_dist1(:) / sum(so4inj_dist1)

    ! Provide diagnostics on the SOA tendencies that affect MXAER.
    do ibin = 1, NBIN
       write(binname, '(A, I2.2)') "MXSOA", ibin

       call addfld(trim(binname)//"CM", (/ 'lev' /), 'A', 'kg/kg/s', 'MXAER SOA gas condensation tendency')
       call addfld(trim(binname)//"PT", (/ 'lev' /), 'A', 'kg/kg/s', 'MXAER SOA photolysis tendency')
    end do

    ! Provide diagnostics for SO4 tendencies from other physics packages
    !
    ! NOTE: This can be useful for determining an SO4 budget and for debugging
    ! SO4 conservation.
    if (carma_do_budget_diags) then

      call addfld("SO4PRBD", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial SO4 pure burden')
      if (carma_diags_file > 0) call add_default("SO4PRBD", carma_diags_file, ' ')
      call addfld("SO4MXBD", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial SO4 mix burden')
      if (carma_diags_file > 0) call add_default("SO4MXBD", carma_diags_file, ' ')
      call addfld("SO4PRCLDBD", horiz_only, 'A', 'kg/m2', 'CARMA, Cloudborne SO4 pure burden')
      if (carma_diags_file > 0) call add_default("SO4PRCLDBD", carma_diags_file, ' ')
      call addfld("SO4MXCLDBD", horiz_only, 'A', 'kg/m2', 'CARMA, Cloudborne SO4 mix burden')

      if (carma_diags_file > 0) call add_default("SO4MXCLDBD", carma_diags_file, ' ')
      call addfld("SO4PRSF", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial SO4 pure surface flux')
      if (carma_diags_file > 0) call add_default("SO4PRSF", carma_diags_file, ' ')
      call addfld("SO4MXSF", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial SO4 mix surface flux')
      if (carma_diags_file > 0) call add_default("SO4MXSF", carma_diags_file, ' ')

      call addfld("H2SO4BD", horiz_only, 'A', 'kg/m2', 'CARMA, H2SO4 burden')
      if (carma_diags_file > 0) call add_default("H2SO4BD", carma_diags_file, ' ')
      call addfld("SO2BD", horiz_only, 'A', 'kg/m2', 'CARMA, SO2 burden')
      if (carma_diags_file > 0) call add_default("SO2BD", carma_diags_file, ' ')

      call addfld("MXBCBD", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial BC mix burden')
      if (carma_diags_file > 0) call add_default("MXBCBD", carma_diags_file, ' ')
      call addfld("MXDUSTBD", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial Dust mix burden')
      if (carma_diags_file > 0) call add_default("MXDUSTBD", carma_diags_file, ' ')
      call addfld("MXOCBD", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial OC mix burden')
      if (carma_diags_file > 0) call add_default("MXOCBD", carma_diags_file, ' ')
      call addfld("MXSALTBD", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial Sea Salt mix burden')
      if (carma_diags_file > 0) call add_default("MXSALTBD", carma_diags_file, ' ')
      call addfld("MXSOABD", horiz_only, 'A', 'kg/m2', 'CARMA, Interstitial SOA mix burden')
      if (carma_diags_file > 0) call add_default("MXSOABD", carma_diags_file, ' ')

      call addfld("MXBCCLDBD", horiz_only, 'A', 'kg/m2', 'CARMA, Cloudborne BC mix burden')
      if (carma_diags_file > 0) call add_default("MXBCCLDBD", carma_diags_file, ' ')
      call addfld("MXDUSTCLDBD", horiz_only, 'A', 'kg/m2', 'CARMA, Cloudborne Dust mix burden')
      if (carma_diags_file > 0) call add_default("MXDUSTCLDBD", carma_diags_file, ' ')
      call addfld("MXOCCLDBD", horiz_only, 'A', 'kg/m2', 'CARMA, Cloudborne OC mix burden')
      if (carma_diags_file > 0) call add_default("MXOCCLDBD", carma_diags_file, ' ')
      call addfld("MXSALTCLDBD", horiz_only, 'A', 'kg/m2', 'CARMA, Cloudborne Sea Salt mix burden')
      if (carma_diags_file > 0) call add_default("MXSALTCLDBD", carma_diags_file, ' ')
      call addfld("MXSOACLDBD", horiz_only, 'A', 'kg/m2', 'CARMA, Cloudborne SOA mix burden')
      if (carma_diags_file > 0) call add_default("MXSOACLDBD", carma_diags_file, ' ')
    end if

    if (carma_do_package_diags) then

      ! Iterate of the packages that have be instrumented. These should match the calls
      ! in physpkg.f90.
      do i = 1, carma_ndiagpkgs
        call addfld("SO4PRBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', SO4 pure burden')
        if (carma_diags_file > 0) call add_default("SO4PRBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SO4MXBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', SO4 mixed burden')
        if (carma_diags_file > 0) call add_default("SO4MXBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("SO4PRSF_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Surface Flux, SO4 pure tendency')
        if (carma_diags_file > 0) call add_default("SO4PRSF_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SO4MXSF_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Surface Flux, SO4 mix tendency')
        if (carma_diags_file > 0) call add_default("SO4MXSF_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("SO4PRTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', SO4 pure tendency')
        if (carma_diags_file > 0) call add_default("SO4PRTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SO4MXTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', SO4 mixed tendency')
        if (carma_diags_file > 0) call add_default("SO4MXTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("SO4PRCLDBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Cloudborne SO4 pure burden')
        if (carma_diags_file > 0) call add_default("SO4PRCLDBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SO4MXCLDBD_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2', trim(carma_diags_packages(i))//', Cloudborne SO4 mixed burden')
        if (carma_diags_file > 0) call add_default("SO4MXCLDBD_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("SO4PRCLDTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Cloudborne SO4 pure tendency')
        if (carma_diags_file > 0) call add_default("SO4PRCLDTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SO4MXCLDTC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', Cloudborne SO4 mixed tendency')
        if (carma_diags_file > 0) call add_default("SO4MXCLDTC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')

        call addfld("H2SO4TC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', H2SO4 total tendency')
        if (carma_diags_file > 0) call add_default("H2SO4TC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
        call addfld("SO2TC_"//trim(carma_diags_packages(i)), horiz_only, 'A', 'kg/m2/s', trim(carma_diags_packages(i))//', SO2 total tendency')
        if (carma_diags_file > 0) call add_default("SO2TC_"//trim(carma_diags_packages(i)), carma_diags_file, ' ')
      end do
    end if

   ! Provide diagnostics for Mass mixing ration summed over the bins
    call addfld("SO4PRMR", (/ 'lev' /), 'A', 'kg/kg', 'SO4 pure mass mixing ratio')
    call addfld("MXSO4MR", (/ 'lev' /), 'A', 'kg/kg', 'SO4 mixed mass mixing ratio')
    call addfld("MXBCMR", (/ 'lev' /), 'A', 'kg/kg', 'BC mixed mass mixing ratio')
    call addfld("MXDUSTMR", (/ 'lev' /), 'A', 'kg/kg', 'DUST mixed mass mixing ratio')
    call addfld("MXOCMR", (/ 'lev' /), 'A', 'kg/kg', 'OC mixed mass mixing ratio')
    call addfld("MXSALTMR", (/ 'lev' /), 'A', 'kg/kg', 'SALT mixed mass mixing ratio')
    call addfld("MXSOAMR", (/ 'lev' /), 'A', 'kg/kg', 'SOA mixed mass mixing ratio')

    return
  end subroutine CARMAMODEL_InitializeModel


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMAMODEL_InitializeParticle(carma, ielem, ibin, latvals, lonvals, mask, q, rc)

    type(carma_type), intent(in)  :: carma      !! the carma object
    integer,          intent(in)  :: ielem      !! element index
    integer,          intent(in)  :: ibin       !! bin index
    real(r8),         intent(in)  :: latvals(:) !! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) !! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    !! Only initialize where .true.
    real(r8),         intent(inout) :: q(:,:)     !! mass mixing ratio (gcol, lev)
    integer,          intent(out) :: rc         !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    !
    ! NOTE: Initialized to 0. by the caller, so nothing needs to be done.

    return
  end subroutine CARMAMODEL_InitializeParticle


  !! This routine is an extension of CARMA_CreateOpticsFile() that allows for
  !! model specific tables to be created in addition to the model independent
  !! methods that are in carma_intr.F90.
  !!
  !! The opticsType that is specified for the group determines how the optical
  !! properties will be generated for that group. Each group can use a different
  !! optics method if needed. Refractive indices need for these calculation are
  !! are specified in the group's elements rather than at the group level. This
  !! allows various mixing approaches to be used to determine the refractive index
  !! for the particle as a whole. If the refractive index for water is needed,
  !! it is specific the the CARMAGAS object for H2O.
  !!
  !! The I_OPTICS_MIXED_YU2105  and I_OPTICS_SULFATE_YU2015 optics methods are
  !! designed to trop_strat models as define in the Yu et al. (2015) paper. The
  !! I_OPTICS_MIXED_YU_H2O includes volume mixing of the water into the shell.
  subroutine CARMAMODEL_CreateOpticsFile(carma, igroup, opticsType, rc)

    implicit none

    type(carma_type), intent(inout)     :: carma         !! the carma object
    integer, intent(in)                 :: igroup        !! group identifier
    integer, intent(in)                 :: opticsType    !! optics type (see I_OPTICS_... in carma_enums.F90)
    integer, intent(out)                :: rc            !! return code, negative indicates failure

    ! Local variables
    logical                             :: do_mie
    integer                             :: cnsttype               ! constituent type

    ! Assume success.
    rc = 0

    ! What type of calculation is needed for this group?
    !
    ! NOTE: Some of these calculations generate optical properties as single mass
    ! coefficients, while others are lookup tables designed around multiple
    ! dimensions.
    select case (opticsType)

      ! This is for the mixed aerosol group as implemented by Yu et al. (2015),
      ! and is specific to the aerosol defintion in that model. There are multiple
      ! elements, some grouped in the core and others in the shell. The refractive
      ! index for the shell is assumed to be only sulfates, and the refractive
      ! index of the core is a mix of dust and black carbon. Core/shell optics
      ! are used to determine the optical properties.
      case(I_OPTICS_MIXED_YU2015)
        call CARMAMODEL_CreateOpticsFile_MixedYu(carma, igroup, rc)
        if (rc < 0) call endrun('carma_CreateOpticsFile::CreateOpticsFile_MixedYu failed.')

      ! This is for the pure sulfate group as implemented by Yu et al. (2015).
      ! The particle may swell, but the refractive index is fixed regardless
      ! of the weight percent of H21SO4 in the particle.
      case(I_OPTICS_SULFATE_YU2015)
        call CARMAMODEL_CreateOpticsFile_SulfateYu(carma, igroup, rc)
        if (rc < 0) call endrun('carma_CreateOpticsFile::CreateOpticsFile_SulfateYu failed.')

      ! This is similar to I_OPTICS_MIXED_YU2015, except that the shell is a volume
      ! mixture of water and H2SO4 rather than just being H2SO4.
      case(I_OPTICS_MIXED_YU_H2O)
        call CARMAMODEL_CreateOpticsFile_MixedYuH2o(carma, igroup, rc)
        if (rc < 0) call endrun('carma_CreateOpticsFile::CreateOpticsFile_MixedYuH2o failed.')

      case default
        call endrun('carma_CreateOpticsFile:: Unknown optics type.')
    end select

    return
  end subroutine CARMAMODEL_CreateOpticsFile


  !! This routine creates files containing optical properties for the mixed group
  !! following Yu et al. (2015). These optical properties are used by the RRTMG radiation
  !! code to include the impact of CARMA particles in the radiative transfer
  !! calculation.
  subroutine CARMAMODEL_CreateOpticsFile_MixedYu(carma, igroup, rc)
    use radconstants, only : nswbands, nlwbands
    use wrap_nf
    use wetr, only         : getwetr

    implicit none

    type(carma_type), intent(inout)     :: carma         !! the carma object
    integer, intent(in)                 :: igroup        !! group index
    integer, intent(out)                :: rc            !! return code, negative indicates failure

    !! Core-shell mixing method for mie and radiation calculations for the Yu et al. (2015)
    !! style table. The CAM optics code will interpolate based upon the current core/shell
    !! mass ratio from a table built using the specified core/shell.
    integer, parameter                  :: ncoreshellratio  = 9               !! Number of core/shell ratio for mie calculations
    integer, parameter                  :: ndstbcratio = 8
    integer, parameter                  :: nkap = 9

    real(kind=f), parameter :: coreshellratio(ncoreshellratio) = (/ 0.001_f, 0.00237_f, 0.00562_f, 0.01333_f, &
                                                                    0.03162_f, 0.07499_f, 0.17782_f, 0.42169_f, 1.0_f /)
    real(kind=f), parameter :: dstbcratio(ndstbcratio) = (/ 0.01_f, 0.025_f, 0.063_f, 0.1_f, 0.3_f, 0.5_f, 0.7_f, 0.9_f/)
    real(kind=f), parameter :: kap(nkap) = (/ 0.1_f, 0.2_f, 0.3_f, 0.4_f, 0.5_f, 0.7_f, 0.9_f, 1.1_f, 1.2_f/)

    ! Local variables
    integer                             :: ibin, iwave, irh, icsr, idb, ikap, icore, ncore
    integer                             :: icorelem(NELEM)
    integer                             :: irhswell
    integer                             :: imiertn
    integer                             :: ienconc
    real(kind=f)                        :: rho(NBIN), rhopwet
    real(kind=f)                        :: r(NBIN), rmass(NBIN), rlow(NBIN), rup(NBIN)
    real(kind=f)                        :: wave(NWAVE)
    complex(kind=f)                     :: refidx(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxS(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxB(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxD(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxC
    !real(kind=f) :: coreimagidx
    character(len=CARMA_NAME_LEN)       :: name
    character(len=CARMA_SHORT_NAME_LEN) :: shortname
    logical                             :: do_mie
    integer                             :: fid
    integer                             :: rhdim, lwdim, swdim, csrdim, dstbcrdim, kapdim
    integer                             :: rhvar, lwvar, swvar, csr_var, dstbcr_var, kap_var
    integer                             :: abs_lw_coreshell_var, qabs_lw_coreshell_var
    integer                             :: ext_sw_coreshell_var, ssa_sw_coreshell_var
    integer                             :: asm_sw_coreshell_var, qext_sw_coreshell_var
    integer                             :: rwetvar
    integer                             :: omdim, andim, namedim
    integer                             :: omvar, anvar, namevar
    integer                             :: dimids(5)
    integer                             :: denvar, slogvar, dryrvar, rminvar, rmaxvar, hygrovar, ntmvar
    real(kind=f)                        :: abs_lw_coreshell(NMIE_RH, nlwbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: qabs_lw_coreshell(NMIE_RH, nlwbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: ext_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: qext_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: ssa_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: asm_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: rwetbin(NMIE_RH)
    character(len=8)                    :: c_name                   ! constituent name
    character(len=32)                   :: aer_name                 ! long enough for both aername and name
    character(len=255)                  :: filepath
    real(kind=f)                        :: rwet
    real(kind=f)                        :: rcore                ! CORE radius used in MIE calculation
    real(kind=f)                        :: Qext
    real(kind=f)                        :: Qsca
    real(kind=f)                        :: asym
    integer                             :: start_text(2), count_text(2)
    integer                             :: sw_r_refidx_var, sw_i_refidx_var, lw_r_refidx_var, lw_i_refidx_var
    integer                             :: ncsr, ndbr
    integer                             :: cnsttype               ! constituent type
    integer                             :: maxbin                 ! last prognostic bin
    integer                             :: LUNOPRT              ! logical unit number for output
    logical                             :: do_print             ! do print output?
    integer                             :: ret

    character(len=32) :: elementname

    ! Assume success.
    rc = 0

    ! Get the wavelength structure.
    call CARMA_GET(carma, rc, wave=wave, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMA_Get failed.')

    ! Get the necessary group properties.
    call CARMAGROUP_Get(carma, igroup, rc, do_mie=do_mie, name=name, shortname=shortname, r=r, &
                        rlow=rlow, rup=rup, rmass=rmass, irhswell=irhswell, imiertn=imiertn, &
                        ienconc=ienconc, ncore=ncore, icorelem=icorelem, cnsttype=cnsttype, maxbin=maxbin)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAGROUP_Get failed.')

    ! The concentration element has the sulfate refractive index.
    call CARMAELEMENT_Get(carma, ienconc, rc, rho=rho, refidx=refidxS)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAELEMENT_Get failed.')

    ! Need to find the dust and black carbon refractive indicies for the core.
    do icore = 1, ncore
      call CARMAELEMENT_Get(carma, icorelem(icore), rc, shortname=elementname, refidx=refidx)
      if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAELEMENT_Get failed.')

      if (trim(elementname) == 'MXBC') then
        refidxB = refidx
      else if (trim(elementname) == 'MXDUST') then
        refidxD = refidx
      end if
    end do


    ! A file needs to be created for each bin.
    do ibin = 1, NBIN

      ! Bins past maxbin are treated as diagnostic even if the group
      ! is prognostic and thus are not advected in the paerent model.
      if (ibin <= maxbin) then

        write(c_name, '(A, I2.2)') trim(shortname), ibin

        ! Construct the path to the file. Each model will have its own subdirectory
        ! where the optical property files are stored.
        filepath = trim(carma_model) // '_' // trim(c_name) // '_rrtmg.nc'

        if (do_print) write(LUNOPRT,*) 'Creating CARMA optics file ... ', trim(filepath)

        ! Create the file.
        call wrap_create(filepath, NF90_CLOBBER, fid)

        ncsr = ncoreshellratio
        ndbr = ndstbcratio

        ! Define the dimensions: rh, lwbands, swbands
        call wrap_def_dim(fid, 'rh_idx',  NMIE_RH,  rhdim)
        call wrap_def_dim(fid, 'lw_band', nlwbands, lwdim)
        call wrap_def_dim(fid, 'sw_band', nswbands, swdim)

        call wrap_def_dim(fid, 'coreshellratio', ncsr, csrdim)
        call wrap_def_dim(fid, 'dstbcratio', ndbr, dstbcrdim)
        call wrap_def_dim(fid, 'kap', nkap, kapdim)

        dimids(1) = rhdim
        call wrap_def_var(fid, 'rh',  NF90_DOUBLE, 1, dimids(1), rhvar)
        call wrap_def_var(fid, 'rwet',NF90_DOUBLE, 1, dimids(1), rwetvar)

        dimids(1) = lwdim
        call wrap_def_var(fid, 'lw_band', NF90_DOUBLE, 1, dimids(1), lwvar)

        dimids(1) = swdim
        call wrap_def_var(fid, 'sw_band', NF90_DOUBLE, 1, dimids(1), swvar)

        dimids(1) = csrdim
        call wrap_def_var(fid, 'coreshellratio', NF90_DOUBLE, 1, dimids(1), csr_var)
        dimids(1) = dstbcrdim
        call wrap_def_var(fid, 'dstbcratio', NF90_DOUBLE, 1, dimids(1), dstbcr_var)
        dimids(1) = kapdim
        call wrap_def_var(fid, 'kap', NF90_DOUBLE, 1, dimids(1), kap_var)


        call wrap_put_att_text(fid, rhvar, 'units', 'fraction')
        call wrap_put_att_text(fid, rwetvar, 'units', 'cm')
        call wrap_put_att_text(fid, lwvar, 'units', 'm')
        call wrap_put_att_text(fid, swvar, 'units', 'm')

        call wrap_put_att_text(fid, csr_var,'units', 'fraction')
        call wrap_put_att_text(fid, dstbcr_var,'units', 'fraction')
        call wrap_put_att_text(fid, kap_var,'units', 'unitless')
        call wrap_put_att_text(fid, csr_var,'long_name', 'coreshell ratio')
        call wrap_put_att_text(fid, dstbcr_var,'long_name', 'dust-bc ratio')
        call wrap_put_att_text(fid, kap_var,'long_name', 'kappa value')

        call wrap_put_att_text(fid, rhvar, 'long_name', 'relative humidity')
        call wrap_put_att_text(fid, rwetvar, 'long_name', 'wet radius')
        call wrap_put_att_text(fid, lwvar, 'long_name', 'longwave bands')
        call wrap_put_att_text(fid, swvar, 'long_name', 'shortwave bands')

        ! Define 3-dimension (:nrh,:nswbands,:ncoreshellratio) LW optics properties: abs_lw_coreshell, qabs_lw_coreshell
        dimids(1) = rhdim
        dimids(2) = lwdim
        dimids(3) = csrdim
        dimids(4) = dstbcrdim
        dimids(5) = kapdim
        call wrap_def_var(fid, 'abs_lw_coreshell', NF90_DOUBLE, 5, dimids(1:5), abs_lw_coreshell_var)
        call wrap_def_var(fid, 'qabs_lw_coreshell',NF90_DOUBLE, 5, dimids(1:5), qabs_lw_coreshell_var)

        call wrap_put_att_text(fid, abs_lw_coreshell_var, 'units', 'meter^2 kilogram^-1')
        call wrap_put_att_text(fid, qabs_lw_coreshell_var,'units', '-')

        ! Define 3-dimension (:nrh,:nswbands,:ncoreshellratio) SW optics properties:
        !  ext_sw_coreshell, qext_sw_coreshell, ssa_sw_coreshell, asm_sw_coreshell
        dimids(1) = rhdim
        dimids(2) = swdim
        dimids(3) = csrdim
        dimids(4) = dstbcrdim
        dimids(5) = kapdim
        call wrap_def_var(fid, 'ext_sw_coreshell', NF90_DOUBLE, 5, dimids(1:5), ext_sw_coreshell_var)
        call wrap_def_var(fid, 'qext_sw_coreshell',NF90_DOUBLE, 5, dimids(1:5), qext_sw_coreshell_var)
        call wrap_def_var(fid, 'ssa_sw_coreshell', NF90_DOUBLE, 5, dimids(1:5), ssa_sw_coreshell_var)
        call wrap_def_var(fid, 'asm_sw_coreshell', NF90_DOUBLE, 5, dimids(1:5), asm_sw_coreshell_var)

        call wrap_put_att_text(fid, ssa_sw_coreshell_var, 'units', 'fraction')
        call wrap_put_att_text(fid, ext_sw_coreshell_var, 'units', 'meter^2 kilogram^-1')
        call wrap_put_att_text(fid, qext_sw_coreshell_var,'units', '-')
        call wrap_put_att_text(fid, asm_sw_coreshell_var, 'units', '-')

        ! Define the variables for the refractive indicies.
        dimids(1) = swdim
        call wrap_def_var(fid, 'refindex_real_aer_sw', NF90_DOUBLE, 1, dimids(1), sw_r_refidx_var)
        call wrap_def_var(fid, 'refindex_im_aer_sw',   NF90_DOUBLE, 1, dimids(1), sw_i_refidx_var)

        dimids(1) = lwdim
        call wrap_def_var(fid, 'refindex_real_aer_lw', NF90_DOUBLE, 1, dimids(1), lw_r_refidx_var)
        call wrap_def_var(fid, 'refindex_im_aer_lw',   NF90_DOUBLE, 1, dimids(1), lw_i_refidx_var)

        call wrap_put_att_text(fid, sw_r_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, sw_i_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, lw_r_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, lw_i_refidx_var, 'units', '-')

        call wrap_put_att_text(fid, sw_r_refidx_var, 'long_name', 'real refractive index of aerosol - shortwave')
        call wrap_put_att_text(fid, sw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - shortwave')
        call wrap_put_att_text(fid, lw_r_refidx_var, 'long_name', 'real refractive index of aerosol - longwave')
        call wrap_put_att_text(fid, lw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - longwave')

        ! Define fields that define the aerosol properties.
        call wrap_def_dim(fid, 'opticsmethod_len',  32, omdim)
        dimids(1) = omdim
        call wrap_def_var(fid, 'opticsmethod',  NF90_CHAR, 1, dimids(1), omvar)

        call wrap_def_dim(fid, 'namelength',  20, andim)
        dimids(1) = andim
        call wrap_def_var(fid, 'aername',  NF90_CHAR, 1, dimids(1), anvar)

        call wrap_def_dim(fid, 'name_len',  32, namedim)
        dimids(1) = namedim
        call wrap_def_var(fid, 'name',  NF90_CHAR, 1, dimids, namevar)

        call wrap_def_var(fid, 'density',            NF90_DOUBLE, 0, dimids(1), denvar)
        call wrap_def_var(fid, 'sigma_logr',         NF90_DOUBLE, 0, dimids(1), slogvar)
        call wrap_def_var(fid, 'dryrad',             NF90_DOUBLE, 0, dimids(1), dryrvar)
        call wrap_def_var(fid, 'radmin_aer',         NF90_DOUBLE, 0, dimids(1), rminvar)
        call wrap_def_var(fid, 'radmax_aer',         NF90_DOUBLE, 0, dimids(1), rmaxvar)
        call wrap_def_var(fid, 'hygroscopicity',     NF90_DOUBLE, 0, dimids(1), hygrovar)
        call wrap_def_var(fid, 'num_to_mass_ratio',  NF90_DOUBLE, 0, dimids(1), ntmvar)

        call wrap_put_att_text(fid, denvar,   'units', 'kg m^-3')
        call wrap_put_att_text(fid, slogvar,  'units', '-')
        call wrap_put_att_text(fid, dryrvar,  'units', 'm')
        call wrap_put_att_text(fid, rminvar,  'units', 'm')
        call wrap_put_att_text(fid, rmaxvar,  'units', 'm')
        call wrap_put_att_text(fid, hygrovar, 'units', '-')
        call wrap_put_att_text(fid, ntmvar,   'units', 'kg^-1')

        call wrap_put_att_text(fid, denvar,   'long_name', 'aerosol material density')
        call wrap_put_att_text(fid, slogvar,  'long_name', 'geometric standard deviation of aerosol')
        call wrap_put_att_text(fid, dryrvar,  'long_name', 'dry number mode radius of aerosol')
        call wrap_put_att_text(fid, rminvar,  'long_name', 'minimum dry radius of aerosol for bin')
        call wrap_put_att_text(fid, rmaxvar,  'long_name', 'maximum dry radius of aerosol for bin')
        call wrap_put_att_text(fid, hygrovar, 'long_name', 'hygroscopicity of aerosol')
        call wrap_put_att_text(fid, ntmvar,   'long_name', 'ratio of number to mass of aerosol')

        ! End the defintion phase of the netcdf file.
        call wrap_enddef(fid)

        ! Write out the dimensions.
        call wrap_put_var_realx(fid, rhvar, mie_rh(:NMIE_RH))
        call wrap_put_var_realx(fid, lwvar, wave(:nlwbands) * 1e-2_f)
        call wrap_put_var_realx(fid, swvar, wave(nlwbands+1:) * 1e-2_f)

        call wrap_put_var_realx(fid, csr_var,coreshellratio(:ncsr))
        call wrap_put_var_realx(fid, dstbcr_var,dstbcratio(:ndstbcratio))
        call wrap_put_var_realx(fid, kap_var,kap(:nkap))

        ! Write out the refractive indicies.
        call wrap_put_var_realx(fid, sw_r_refidx_var, real(refidxS(nlwbands+1:, 1)))
        call wrap_put_var_realx(fid, sw_i_refidx_var, aimag(refidxS(nlwbands+1:, 1)))
        call wrap_put_var_realx(fid, lw_r_refidx_var, real(refidxS(:nlwbands, 1)))
        call wrap_put_var_realx(fid, lw_i_refidx_var, aimag(refidxS(:nlwbands, 1)))

        ! Pad the names out with spaces.
        aer_name = '                                '
        aer_name(1:len(trim(c_name))) = c_name

        start_text(1) = 1
        count_text(1) = 32
        call wrap_put_vara_text(fid, namevar, start_text, count_text, (/ aer_name /))
        count_text(1) = 20
        call wrap_put_vara_text(fid, anvar, start_text, count_text, (/ aer_name /))

        count_text(1) = len('hygroscopic_coreshell           ')
        call wrap_put_vara_text(fid, omvar, start_text, count_text, (/ 'hygroscopic_coreshell           ' /))

        call wrap_put_var_realx(fid, denvar,   (/ rho(ibin) * 1e-3_f / 1e-6_f /))
        call wrap_put_var_realx(fid, slogvar,  (/ 0._f /))
        call wrap_put_var_realx(fid, dryrvar,  (/ r(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, rminvar,  (/ rlow(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, rmaxvar,  (/ rup(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, hygrovar, (/ 0.6_f /))
        call wrap_put_var_realx(fid, ntmvar,   (/ 1._f / rmass(ibin) / 1e-3_f /))

        ! For now, ext_sw(:nrh, :nswbands) and ext_sw_coreshell(:nrh, :nswbands, :ncoreshellratio) both are calculated
        ! Since other aerosols in CAM may use ext_sw rather than ext_sw_coreshell
        ! Modified by Pengfei Yu
        ! April.1, 2012

        !--------------------------- for 5-D core-shell optical properties ----------------------------

        ! Iterate over a range of relative humidities, since the particle may swell
        ! with relative humidity which will change its optical properties.
        do irh = 1, NMIE_RH

          do ikap = 1, nkap

            ! Determine the wet radius.
            call getwetr(carma, igroup, mie_rh(irh), r(ibin), rwet, rho(ibin), rhopwet, rc, kappa=kap(ikap), temp=270._f)
            rwetbin(irh) = rwet

            ! Calculate at each wavelength.
            do iwave = 1, NWAVE

              ! For now just assume BC/OC constant 15%
              ! rcore = r(ibin)*(0.15**onethird)
              ! Using Mie code, consider core/shell ratio
              do icsr = 1, ncsr
                if (ncsr > 1) then
                  rcore = r(ibin)*(coreshellratio(icsr)**onethird)
                else
                  rcore = 0.0_f
                endif

                ! Using Mie code, assume the particle is CORE-SHELL
                ! By: Pengfei Yu
                ! Mar.22, 2012

                !write(*,*) 'before call mie-3D, icsr = ', icsr, ' ;iwave = ', iwave, ' ;irh = ', irh
                !write(*,*) 'ibin = ', ibin, ' ;rcore = ', rcore, ' ;csratio = ', coreshellratio(icsr)

                do idb = 1, ndbr

                  ! NOTE: This is not the best way to combine the dust and BC refractive indices
                  ! for the core. Volume mixing should be used for both the real and imaginary
                  ! parts, not just the imaginary.
!                  coreimagidx = dstbcratio(idb) * aimag(refidxB(iwave,1)) + (1._f - dstbcratio(idb)) * aimag(refidxD(iwave,1))
!                  refidxC = cmplx((real(refidxD(iwave,1)) + real(refidxB(iwave,1))) / 2._f, coreimagidx)
                  refidxC = dstbcratio(idb) * refidxB(iwave,1) + (1._f - dstbcratio(idb)) * refidxD(iwave,1)

                  call mie(carma, &
                           imiertn, &
                           rwet, &
                           wave(iwave), &
                           0._f, &
                           3.0_f, &
                           0.0_f, &
                           1.0_f, &
                           refidxS(iwave, 1), &
                           rcore, &
                           refidxC, &
                           Qext, &
                           Qsca, &
                           asym, &
                           rc)
                  if (rc < 0) call endrun('carma_CreateOpticsFile::mie failed.')

                  ! Calculate  the shortwave and longwave properties?
                  !
                  ! NOTE: miess is in cgs units, but the optics file needs to be in mks
                  ! units, so perform the necessary conversions.
                  if (iwave <= nlwbands) then

                    ! Longwave just needs absorption: abs_lw.
                    qabs_lw_coreshell(irh, iwave, icsr, idb, ikap) = (Qext - Qsca) ! absorption per particle
                    abs_lw_coreshell (irh, iwave, icsr, idb, ikap) = (Qext - Qsca) * PI * (rwet * 1e-2_f)**2 &
                                                                      / (rmass(ibin) * 1e-3_f)
                  else

                    ! Shortwave needs extinction, single scattering albedo and asymmetry factor:
                    ! ext_sw, qext_sw, ssa_sw and asm_sw.
                    qext_sw_coreshell(irh, iwave - nlwbands, icsr, idb, ikap) = Qext ! extinction per particle
                    ext_sw_coreshell (irh, iwave - nlwbands, icsr, idb, ikap) = Qext * PI * (rwet * 1e-2_f)**2 &
                                                                                / (rmass(ibin) * 1e-3_f)
                    ssa_sw_coreshell (irh, iwave - nlwbands, icsr, idb, ikap) = Qsca / Qext
                    asm_sw_coreshell (irh, iwave - nlwbands, icsr, idb, ikap) = asym
                  end if
                end do   ! idb
              end do       ! icsr
            end do   ! iwave
          end do      ! ikap
        end do     ! irh

        call wrap_put_var_realx(fid, rwetvar, rwetbin(:))

        ! Write out the longwave fields.
        ret = nf90_put_var(fid, abs_lw_coreshell_var, abs_lw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', abs_lw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, qabs_lw_coreshell_var,  qabs_lw_coreshell(:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', qabs_lw_coreshell_var
           call handle_error(ret)
        end if

        ! Write out the shortwave fields.
        ret = nf90_put_var(fid, ext_sw_coreshell_var,   ext_sw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', ext_sw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, qext_sw_coreshell_var,  qext_sw_coreshell(:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', qext_sw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, ssa_sw_coreshell_var,   ssa_sw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', ssa_sw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, asm_sw_coreshell_var,   asm_sw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', asm_sw_coreshell_var
           call handle_error(ret)
        end if

        ! Close the file.
        call wrap_close(fid)
      end if
    end do

    return
  end subroutine CARMAMODEL_CreateOpticsFile_MixedYu

  !! This routine creates files containing optical properties for the mixed group
  !! following Yu et al. (2015), except that it includes water vapor in the shell.
  !! The difference between the wet and dry radius is assumed to be water valor and
  !! the shell is a volume mix of the H2SO4 and the water. These optical properties
  !! are used by the RRTMG radiation code to include the impact of CARMA particles
  !! in the radiative transfer calculation.
  !!
  !! NOTE: The table structure is the same as for MixedYu, so no changes need to be
  !! made on the CAM side to use these optics.
  subroutine CARMAMODEL_CreateOpticsFile_MixedYuH2o(carma, igroup, rc)
    use radconstants, only : nswbands, nlwbands
    use wrap_nf
    use wetr, only         : getwetr

    implicit none

    type(carma_type), intent(inout)     :: carma         !! the carma object
    integer, intent(in)                 :: igroup        !! group index
    integer, intent(out)                :: rc            !! return code, negative indicates failure

    !! Core-shell mixing method for mie and radiation calculations for the Yu et al. (2015)
    !! style table. The CAM optics code will interpolate based upon the current core/shell
    !! mass ratio from a table built using the specified core/shell.
    integer, parameter                  :: ncoreshellratio  = 9               !! Number of core/shell ratio for mie calculations
    integer, parameter	                :: ndstbcratio = 8
    integer, parameter	                :: nkap = 9

    real(kind=f)                        :: coreshellratio(ncoreshellratio) = (/ 0.001_f, 0.00237_f, 0.00562_f, 0.01333_f, 0.03162_f, 0.07499_f, 0.17782_f, 0.42169_f, 1.0_f /)
    real(kind=f)		                    :: dstbcratio(ndstbcratio) = (/ 0.01_f, 0.025_f, 0.063_f, 0.1_f, 0.3_f, 0.5_f, 0.7_f, 0.9_f/)
    real(kind=f)		                    :: kap(nkap) = (/ 0.1_f, 0.2_f, 0.3_f, 0.4_f, 0.5_f, 0.7_f, 0.9_f, 1.1_f, 1.2_f/)

    ! Local variables
    integer                             :: ibin, iwave, irh, icsr, idb, ikap, icore, ncore
    integer                             :: icorelem(NELEM)
    integer                             :: irhswell
    integer                             :: imiertn
    integer                             :: ienconc
    real(kind=f)                        :: rho(NBIN), rhopwet
    real(kind=f)                        :: r(NBIN), rmass(NBIN), rlow(NBIN), rup(NBIN)
    real(kind=f)                        :: wave(NWAVE)
    complex(kind=f)                     :: refidx(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxS(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxB(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxD(NWAVE, NREFIDX)
    complex(kind=f)                     :: refidxW(NWAVE)
    complex(kind=f)                     :: refidxC
    complex(kind=f)                     :: refidxSH
    !real(kind=f) :: coreimagidx
    character(len=CARMA_NAME_LEN)       :: name
    character(len=CARMA_SHORT_NAME_LEN) :: shortname
    logical                             :: do_mie
    integer                             :: fid
    integer                             :: rhdim, lwdim, swdim, csrdim, dstbcrdim, kapdim
    integer                             :: rhvar, lwvar, swvar, csr_var, dstbcr_var, kap_var
    integer                             :: abs_lw_coreshell_var, qabs_lw_coreshell_var
    integer                             :: ext_sw_coreshell_var, ssa_sw_coreshell_var, asm_sw_coreshell_var, qext_sw_coreshell_var
    integer                             :: rwetvar
    integer                             :: omdim, andim, namedim
    integer                             :: omvar, anvar, namevar
    integer                             :: dimids(5)
    integer                             :: denvar, slogvar, dryrvar, rminvar, rmaxvar, hygrovar, ntmvar
    real(kind=f)                        :: abs_lw_coreshell(NMIE_RH, nlwbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: qabs_lw_coreshell(NMIE_RH, nlwbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: ext_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: qext_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: ssa_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: asm_sw_coreshell(NMIE_RH, nswbands, ncoreshellratio, ndstbcratio, nkap)
    real(kind=f)                        :: rwetbin(NMIE_RH)
    character(len=8)                    :: c_name                   ! constituent name
    character(len=32)                   :: aer_name                 ! long enough for both aername and name
    character(len=255)                  :: filepath
    real(kind=f)                        :: rwet
    real(kind=f)                        :: rcore		! CORE radius used in MIE calculation
    real(kind=f)                        :: Qext
    real(kind=f)                        :: Qsca
    real(kind=f)                        :: asym
    integer                             :: start_text(2), count_text(2)
    integer                             :: sw_r_refidx_var, sw_i_refidx_var, lw_r_refidx_var, lw_i_refidx_var
    integer                             :: ncsr, ndbr
    integer                             :: cnsttype               ! constituent type
    integer                             :: maxbin                 ! last prognostic bin
    integer                             :: LUNOPRT              ! logical unit number for output
    logical                             :: do_print             ! do print output?
    integer                             :: ret
    real(kind=f)                        :: volwater
    real(kind=f)                        :: volsulfate
    real(kind=f)                        :: volshell
    integer                             :: igash2o

    character(len=32) :: elementname

    ! Assume success.
    rc = 0

    ! Get the wavelength structure.
    call CARMA_GET(carma, rc, wave=wave, do_print=do_print, LUNOPRT=LUNOPRT, igash2o=igash2o)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMA_Get failed.')

    ! Get the necessary group properties.
    call CARMAGROUP_Get(carma, igroup, rc, do_mie=do_mie, name=name, shortname=shortname, r=r, &
                        rlow=rlow, rup=rup, rmass=rmass, irhswell=irhswell, imiertn=imiertn, &
                        ienconc=ienconc, ncore=ncore, icorelem=icorelem, cnsttype=cnsttype, maxbin=maxbin)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAGROUP_Get failed.')

    ! The concentration element has the sulfate refractive index.
    call CARMAELEMENT_Get(carma, ienconc, rc, rho=rho, refidx=refidxS)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAELEMENT_Get failed.')

    ! Need to find the dust and black carbon refractive indicies for the core.
    do icore = 1, ncore
      call CARMAELEMENT_Get(carma, icorelem(icore), rc, shortname=elementname, refidx=refidx)
      if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAELEMENT_Get failed.')

      if (trim(elementname) == 'MXBC') then
        refidxB = refidx
      else if (trim(elementname) == 'MXDUST') then
        refidxD = refidx
      end if
    end do

    ! Get the refractive index for water.
    call CARMAGAS_Get(carma, igash2o, rc, refidx=refidxW)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAGAS_Get failed.')

    refidxW(:) = CMPLX(waterreal(:), waterimag(:), kind=f)

    ! A file needs to be created for each bin.
    do ibin = 1, NBIN

      ! Bins past maxbin are treated as diagnostic even if the group
      ! is prognostic and thus are not advected in the paerent model.
      if (ibin <= maxbin) then

        write(c_name, '(A, I2.2)') trim(shortname), ibin

        ! Construct the path to the file. Each model will have its own subdirectory
        ! where the optical property files are stored.
        filepath = trim(carma_model) // '_' // trim(c_name) // '_rrtmg.nc'

        if (do_print) write(LUNOPRT,*) 'Creating CARMA optics file ... ', trim(filepath)

        ! Create the file.
        call wrap_create(filepath, NF90_CLOBBER, fid)

        ncsr = ncoreshellratio
        ndbr = ndstbcratio

        ! Define the dimensions: rh, lwbands, swbands
        call wrap_def_dim(fid, 'rh_idx',  NMIE_RH,  rhdim)
        call wrap_def_dim(fid, 'lw_band', nlwbands, lwdim)
        call wrap_def_dim(fid, 'sw_band', nswbands, swdim)

        call wrap_def_dim(fid, 'coreshellratio', ncsr, csrdim)
        call wrap_def_dim(fid, 'dstbcratio', ndbr, dstbcrdim)
        call wrap_def_dim(fid, 'kap', nkap, kapdim)

        dimids(1) = rhdim
        call wrap_def_var(fid, 'rh',  NF90_DOUBLE, 1, dimids(1), rhvar)
        call wrap_def_var(fid, 'rwet',NF90_DOUBLE, 1, dimids(1), rwetvar)

        dimids(1) = lwdim
        call wrap_def_var(fid, 'lw_band', NF90_DOUBLE, 1, dimids(1), lwvar)

        dimids(1) = swdim
        call wrap_def_var(fid, 'sw_band', NF90_DOUBLE, 1, dimids(1), swvar)

        dimids(1) = csrdim
        call wrap_def_var(fid, 'coreshellratio', NF90_DOUBLE, 1, dimids(1), csr_var)
        dimids(1) = dstbcrdim
        call wrap_def_var(fid, 'dstbcratio', NF90_DOUBLE, 1, dimids(1), dstbcr_var)
        dimids(1) = kapdim
        call wrap_def_var(fid, 'kap', NF90_DOUBLE, 1, dimids(1), kap_var)


        call wrap_put_att_text(fid, rhvar, 'units', 'fraction')
        call wrap_put_att_text(fid, rwetvar, 'units', 'cm')
        call wrap_put_att_text(fid, lwvar, 'units', 'm')
        call wrap_put_att_text(fid, swvar, 'units', 'm')

        call wrap_put_att_text(fid, csr_var,'units', 'fraction')
        call wrap_put_att_text(fid, dstbcr_var,'units', 'fraction')
        call wrap_put_att_text(fid, kap_var,'units', 'unitless')
        call wrap_put_att_text(fid, csr_var,'long_name', 'coreshell ratio')
        call wrap_put_att_text(fid, dstbcr_var,'long_name', 'dust-bc ratio')
        call wrap_put_att_text(fid, kap_var,'long_name', 'kappa value')

        call wrap_put_att_text(fid, rhvar, 'long_name', 'relative humidity')
        call wrap_put_att_text(fid, rwetvar, 'long_name', 'wet radius')
        call wrap_put_att_text(fid, lwvar, 'long_name', 'longwave bands')
        call wrap_put_att_text(fid, swvar, 'long_name', 'shortwave bands')

        ! Define 3-dimension (:nrh,:nswbands,:ncoreshellratio) LW optics properties: abs_lw_coreshell, qabs_lw_coreshell
        dimids(1) = rhdim
        dimids(2) = lwdim
        dimids(3) = csrdim
        dimids(4) = dstbcrdim
        dimids(5) = kapdim
        call wrap_def_var(fid, 'abs_lw_coreshell', NF90_DOUBLE, 5, dimids(1:5), abs_lw_coreshell_var)
        call wrap_def_var(fid, 'qabs_lw_coreshell',NF90_DOUBLE, 5, dimids(1:5), qabs_lw_coreshell_var)

        call wrap_put_att_text(fid, abs_lw_coreshell_var, 'units', 'meter^2 kilogram^-1')
        call wrap_put_att_text(fid, qabs_lw_coreshell_var,'units', '-')

        ! Define 3-dimension (:nrh,:nswbands,:ncoreshellratio) SW optics properties:
        !    ext_sw_coreshell, qext_sw_coreshell, ssa_sw_coreshell, asm_sw_coreshell
        dimids(1) = rhdim
        dimids(2) = swdim
        dimids(3) = csrdim
        dimids(4) = dstbcrdim
        dimids(5) = kapdim
        call wrap_def_var(fid, 'ext_sw_coreshell', NF90_DOUBLE, 5, dimids(1:5), ext_sw_coreshell_var)
        call wrap_def_var(fid, 'qext_sw_coreshell',NF90_DOUBLE, 5, dimids(1:5), qext_sw_coreshell_var)
        call wrap_def_var(fid, 'ssa_sw_coreshell', NF90_DOUBLE, 5, dimids(1:5), ssa_sw_coreshell_var)
        call wrap_def_var(fid, 'asm_sw_coreshell', NF90_DOUBLE, 5, dimids(1:5), asm_sw_coreshell_var)

        call wrap_put_att_text(fid, ssa_sw_coreshell_var, 'units', 'fraction')
        call wrap_put_att_text(fid, ext_sw_coreshell_var, 'units', 'meter^2 kilogram^-1')
        call wrap_put_att_text(fid, qext_sw_coreshell_var,'units', '-')
        call wrap_put_att_text(fid, asm_sw_coreshell_var, 'units', '-')

        ! Define the variables for the refractive indicies.
        dimids(1) = swdim
        call wrap_def_var(fid, 'refindex_real_aer_sw', NF90_DOUBLE, 1, dimids(1), sw_r_refidx_var)
        call wrap_def_var(fid, 'refindex_im_aer_sw',   NF90_DOUBLE, 1, dimids(1), sw_i_refidx_var)

        dimids(1) = lwdim
        call wrap_def_var(fid, 'refindex_real_aer_lw', NF90_DOUBLE, 1, dimids(1), lw_r_refidx_var)
        call wrap_def_var(fid, 'refindex_im_aer_lw',   NF90_DOUBLE, 1, dimids(1), lw_i_refidx_var)

        call wrap_put_att_text(fid, sw_r_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, sw_i_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, lw_r_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, lw_i_refidx_var, 'units', '-')

        call wrap_put_att_text(fid, sw_r_refidx_var, 'long_name', 'real refractive index of aerosol - shortwave')
        call wrap_put_att_text(fid, sw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - shortwave')
        call wrap_put_att_text(fid, lw_r_refidx_var, 'long_name', 'real refractive index of aerosol - longwave')
        call wrap_put_att_text(fid, lw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - longwave')

        ! Define fields that define the aerosol properties.
        call wrap_def_dim(fid, 'opticsmethod_len',  32, omdim)
        dimids(1) = omdim
        call wrap_def_var(fid, 'opticsmethod',  NF90_CHAR, 1, dimids(1), omvar)

        call wrap_def_dim(fid, 'namelength',  20, andim)
        dimids(1) = andim
        call wrap_def_var(fid, 'aername',  NF90_CHAR, 1, dimids(1), anvar)

        call wrap_def_dim(fid, 'name_len',  32, namedim)
        dimids(1) = namedim
        call wrap_def_var(fid, 'name',  NF90_CHAR, 1, dimids, namevar)

        call wrap_def_var(fid, 'density',            NF90_DOUBLE, 0, dimids(1), denvar)
        call wrap_def_var(fid, 'sigma_logr',         NF90_DOUBLE, 0, dimids(1), slogvar)
        call wrap_def_var(fid, 'dryrad',             NF90_DOUBLE, 0, dimids(1), dryrvar)
        call wrap_def_var(fid, 'radmin_aer',         NF90_DOUBLE, 0, dimids(1), rminvar)
        call wrap_def_var(fid, 'radmax_aer',         NF90_DOUBLE, 0, dimids(1), rmaxvar)
        call wrap_def_var(fid, 'hygroscopicity',     NF90_DOUBLE, 0, dimids(1), hygrovar)
        call wrap_def_var(fid, 'num_to_mass_ratio',  NF90_DOUBLE, 0, dimids(1), ntmvar)

        call wrap_put_att_text(fid, denvar,   'units', 'kg m^-3')
        call wrap_put_att_text(fid, slogvar,  'units', '-')
        call wrap_put_att_text(fid, dryrvar,  'units', 'm')
        call wrap_put_att_text(fid, rminvar,  'units', 'm')
        call wrap_put_att_text(fid, rmaxvar,  'units', 'm')
        call wrap_put_att_text(fid, hygrovar, 'units', '-')
        call wrap_put_att_text(fid, ntmvar,   'units', 'kg^-1')

        call wrap_put_att_text(fid, denvar,   'long_name', 'aerosol material density')
        call wrap_put_att_text(fid, slogvar,  'long_name', 'geometric standard deviation of aerosol')
        call wrap_put_att_text(fid, dryrvar,  'long_name', 'dry number mode radius of aerosol')
        call wrap_put_att_text(fid, rminvar,  'long_name', 'minimum dry radius of aerosol for bin')
        call wrap_put_att_text(fid, rmaxvar,  'long_name', 'maximum dry radius of aerosol for bin')
        call wrap_put_att_text(fid, hygrovar, 'long_name', 'hygroscopicity of aerosol')
        call wrap_put_att_text(fid, ntmvar,   'long_name', 'ratio of number to mass of aerosol')

        ! End the defintion phase of the netcdf file.
        call wrap_enddef(fid)

        ! Write out the dimensions.
        call wrap_put_var_realx(fid, rhvar, mie_rh(:NMIE_RH))
        call wrap_put_var_realx(fid, lwvar, wave(:nlwbands) * 1e-2_f)
        call wrap_put_var_realx(fid, swvar, wave(nlwbands+1:) * 1e-2_f)

        call wrap_put_var_realx(fid, csr_var,coreshellratio(:ncsr))
        call wrap_put_var_realx(fid, dstbcr_var,dstbcratio(:ndstbcratio))
        call wrap_put_var_realx(fid, kap_var,kap(:nkap))

        ! Write out the refractive indicies.
        call wrap_put_var_realx(fid, sw_r_refidx_var, real(refidxS(nlwbands+1:, 1)))
        call wrap_put_var_realx(fid, sw_i_refidx_var, aimag(refidxS(nlwbands+1:, 1)))
        call wrap_put_var_realx(fid, lw_r_refidx_var, real(refidxS(:nlwbands, 1)))
        call wrap_put_var_realx(fid, lw_i_refidx_var, aimag(refidxS(:nlwbands, 1)))

        ! Pad the names out with spaces.
        aer_name = '                                '
        aer_name(1:len(trim(c_name))) = c_name

        start_text(1) = 1
        count_text(1) = 32
        call wrap_put_vara_text(fid, namevar, start_text, count_text, (/ aer_name /))
        count_text(1) = 20
        call wrap_put_vara_text(fid, anvar, start_text, count_text, (/ aer_name /))

        count_text(1) = len('hygroscopic_coreshell           ')
        call wrap_put_vara_text(fid, omvar, start_text, count_text, (/ 'hygroscopic_coreshell           ' /))

        call wrap_put_var_realx(fid, denvar,   (/ rho(ibin) * 1e-3_f / 1e-6_f /))
        call wrap_put_var_realx(fid, slogvar,  (/ 0._f /))
        call wrap_put_var_realx(fid, dryrvar,  (/ r(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, rminvar,  (/ rlow(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, rmaxvar,  (/ rup(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, hygrovar, (/ 0.6_f /))
        call wrap_put_var_realx(fid, ntmvar,   (/ 1._f / rmass(ibin) / 1e-3_f /))

        ! For now, ext_sw(:nrh, :nswbands) and ext_sw_coreshell(:nrh, :nswbands, :ncoreshellratio) both are calculated
        ! Since other aerosols in CAM may use ext_sw rather than ext_sw_coreshell
        ! Modified by Pengfei Yu
        ! April.1, 2012

        !--------------------------- for 5-D core-shell optical properties ----------------------------

        ! Iterate over a range of relative humidities, since the particle may swell
        ! with relative humidity which will change its optical properties.
        do irh = 1, NMIE_RH

          do ikap = 1, nkap

            ! Determine the wet radius.
            call getwetr(carma, igroup, mie_rh(irh), r(ibin), rwet, rho(ibin), rhopwet, rc, kappa=kap(ikap), temp=270._f)
            rwetbin(irh) = rwet

            ! Calculate at each wavelength.
            do iwave = 1, NWAVE

              ! For now just assume BC/OC constant 15%
              ! rcore = r(ibin)*(0.15**onethird)
              ! Using Mie code, consider core/shell ratio
              do icsr = 1, ncsr
                if (ncsr > 1) then
                  rcore = r(ibin)*(coreshellratio(icsr)**onethird)
                else
                  rcore = 0.0_f
                endif

                ! This is not in Yu (2015), but rather than using the refractive
                ! index of H2SO4 for the shell, do a volume mix of water and H2SO4
                ! for the refractive index of the shell.
                volwater = rwet**3._f - r(ibin)**3._f
                volsulfate = r(ibin)**3._f * (1._f - coreshellratio(icsr))
                volshell = volwater + volsulfate
                if (volshell > 0._f) then
                  refidxSH = (volwater / volshell) * refidxW(iwave) + (volsulfate / volshell) * refidxS(iwave, 1)
                else
                  refidxSH = refidxS(iwave, 1)
                end if

                ! Using Mie code, assume the particle is CORE-SHELL
                ! By: Pengfei Yu
                ! Mar.22, 2012

                !write(*,*) 'before call mie-3D, icsr = ', icsr, ' ;iwave = ', iwave, ' ;irh = ', irh
                !write(*,*) 'ibin = ', ibin, ' ;rcore = ', rcore, ' ;csratio = ', coreshellratio(icsr)

                do idb = 1, ndbr

                  ! NOTE: This is not the best way to combine the dust and BC refractive indices
                  ! for the core. Volume mixing should be used for both the real and imaginary
                  ! parts, not just the imaginary.
!                  coreimagidx = dstbcratio(idb) * aimag(refidxB(iwave,1)) + (1._f - dstbcratio(idb)) * aimag(refidxD(iwave,1))
!                  refidxC = cmplx((real(refidxD(iwave,1)) + real(refidxB(iwave,1))) / 2._f, coreimagidx)
                  refidxC = dstbcratio(idb) * refidxB(iwave,1) + (1._f - dstbcratio(idb)) * refidxD(iwave,1)

                  call mie(carma, &
                           imiertn, &
                           rwet, &
                           wave(iwave), &
                           0._f, &
                           3.0_f, &
                           0.0_f, &
                           1.0_f, &
                           refidxSH, &
                           rcore, &
                           refidxC, &
                           Qext, &
                           Qsca, &
                           asym, &
                           rc)
                  if (rc < 0) call endrun('carma_CreateOpticsFile::mie failed.')

                  ! Calculate  the shortwave and longwave properties?
                  !
                  ! NOTE: miess is in cgs units, but the optics file needs to be in mks
                  ! units, so perform the necessary conversions.
                  if (iwave <= nlwbands) then

                    ! Longwave just needs absorption: abs_lw.
                    qabs_lw_coreshell(irh, iwave, icsr, idb, ikap) = (Qext - Qsca)                            ! absorption per particle
                    abs_lw_coreshell (irh, iwave, icsr, idb, ikap) = (Qext - Qsca) * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f)
                  else

                    ! Shortwave needs extinction, single scattering albedo and asymmetry factor:
                    ! ext_sw, qext_sw, ssa_sw and asm_sw.
                    qext_sw_coreshell(irh, iwave - nlwbands, icsr, idb, ikap) = Qext                          ! extinction per particle
                    ext_sw_coreshell (irh, iwave - nlwbands, icsr, idb, ikap) = Qext * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f)
                    ssa_sw_coreshell (irh, iwave - nlwbands, icsr, idb, ikap) = Qsca / Qext
                    asm_sw_coreshell (irh, iwave - nlwbands, icsr, idb, ikap) = asym
                  end if
                end do   ! idb
              end do       ! icsr
            end do   ! iwave
          end do      ! ikap
        end do     ! irh

        call wrap_put_var_realx(fid, rwetvar, rwetbin(:))

        ! Write out the longwave fields.
        ret = nf90_put_var(fid, abs_lw_coreshell_var, abs_lw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', abs_lw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, qabs_lw_coreshell_var,  qabs_lw_coreshell(:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', qabs_lw_coreshell_var
           call handle_error(ret)
        end if

        ! Write out the shortwave fields.
        ret = nf90_put_var(fid, ext_sw_coreshell_var,   ext_sw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', ext_sw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, qext_sw_coreshell_var,  qext_sw_coreshell(:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', qext_sw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, ssa_sw_coreshell_var,   ssa_sw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', ssa_sw_coreshell_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, asm_sw_coreshell_var,   asm_sw_coreshell (:, :, :, :, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_MixedYu: error writing varid =', asm_sw_coreshell_var
           call handle_error(ret)
        end if

        ! Close the file.
        call wrap_close(fid)
      end if
    end do

    return
  end subroutine CARMAMODEL_CreateOpticsFile_MixedYuH2o


  !! This routine creates files containing optical properties for the pure sulfate group
  !! following Yu et al. (2015). These optical properties are used by the RRTMG radiation
  !! code to include the impact of CARMA particles in the radiative transfer
  !! calculation.
  subroutine CARMAMODEL_CreateOpticsFile_SulfateYu(carma, igroup, rc)
    use radconstants, only : nswbands, nlwbands
    use wrap_nf
    use wetr, only         : getwetr

    implicit none

    type(carma_type), intent(inout)     :: carma         !! the carma object
    integer, intent(in)                 :: igroup        !! group index
    integer, intent(out)                :: rc            !! return code, negative indicates failure

    ! Local variables
    integer                             :: ibin, iwave, iwtp
    integer                             :: irhswell
    integer                             :: imiertn
    integer                             :: ienconc
    real(kind=f)                        :: rho(NBIN), rhopwet
    real(kind=f)                        :: r(NBIN), rmass(NBIN), rlow(NBIN), rup(NBIN)
    real(kind=f)                        :: wave(NWAVE)
    complex(kind=f)                     :: refidx(NWAVE, NREFIDX)
    character(len=CARMA_NAME_LEN)       :: name
    character(len=CARMA_SHORT_NAME_LEN) :: shortname
    integer                             :: fid
    integer                             :: rhdim, lwdim, swdim, wtpdim
    integer                             :: rhvar, lwvar, swvar, wtp_var
    integer                             :: rwetvar
    integer				:: abs_lw_wtp_var, qabs_lw_wtp_var
    integer                             :: ext_sw_wtp_var, ssa_sw_wtp_var, asm_sw_wtp_var, qext_sw_wtp_var
    integer                             :: omdim, andim, namedim
    integer                             :: omvar, anvar, namevar
    integer                             :: dimids(2)
    integer                             :: denvar, slogvar, dryrvar, rminvar, rmaxvar, hygrovar, ntmvar
    real(kind=f)                        :: abs_lw_wtp(NMIE_WTP, nlwbands)
    real(kind=f)                        :: qabs_lw_wtp(NMIE_WTP, nlwbands)
    real(kind=f)                        :: ext_sw_wtp(NMIE_WTP, nswbands)
    real(kind=f)                        :: qext_sw_wtp(NMIE_WTP, nswbands)
    real(kind=f)                        :: ssa_sw_wtp(NMIE_WTP, nswbands)
    real(kind=f)                        :: asm_sw_wtp(NMIE_WTP, nswbands)
    character(len=8)                    :: c_name                   ! constituent name
    character(len=32)                   :: aer_name                 ! long enough for both aername and name
    character(len=255)                  :: filepath
    real(kind=f)                        :: rwet
    real(kind=f)                        :: Qext
    real(kind=f)                        :: Qsca
    real(kind=f)                        :: asym
    integer                             :: start_text(2), count_text(2)
    integer                             :: sw_r_refidx_var, sw_i_refidx_var, lw_r_refidx_var, lw_i_refidx_var
    integer                             :: cnsttype               ! constituent type
    integer                             :: maxbin                 ! last prognostic bin
    integer                             :: LUNOPRT              ! logical unit number for output
    logical                             :: do_print             ! do print output?
    integer                             :: ret


    ! Assume success.
    rc = 0

    ! Get the wavelength structure.
    call CARMA_GET(carma, rc, wave=wave, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMA_Get failed.')

    ! Get the necessary group properties.
    call CARMAGROUP_Get(carma, igroup, rc, name=name, shortname=shortname, r=r, &
                        rlow=rlow, rup=rup, rmass=rmass, irhswell=irhswell, &
                        ienconc=ienconc, cnsttype=cnsttype, maxbin=maxbin, imiertn=imiertn)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAGROUP_Get failed.')

    ! Get the necessary element properties.
    call CARMAELEMENT_Get(carma, ienconc, rc, rho=rho, refidx=refidx)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAELEMENT_Get failed.')

    ! A file needs to be created for each bin.
    do ibin = 1, NBIN

      ! Bins past maxbin are treated as diagnostic even if the group
      ! is prognostic and thus are not advected in the paerent model.
      if (ibin <= maxbin) then

        write(c_name, '(A, I2.2)') trim(shortname), ibin

        ! Construct the path to the file. Each model will have its own subdirectory
        ! where the optical property files are stored.
        filepath = trim(carma_model) // '_' // trim(c_name) // '_rrtmg.nc'

        if (do_print) write(LUNOPRT,*) 'Creating CARMA optics file ... ', trim(filepath)

        ! Create the file.
        call wrap_create(filepath, NF90_CLOBBER, fid)

        ! Define the dimensions: rh, lwbands, swbands
        call wrap_def_dim(fid, 'rh_idx',  NMIE_RH,  rhdim)
        call wrap_def_dim(fid, 'lw_band', nlwbands, lwdim)
        call wrap_def_dim(fid, 'sw_band', nswbands, swdim)

        call wrap_def_dim(fid, 'wgtpct', NMIE_WTP, wtpdim)

        dimids(1) = rhdim
        call wrap_def_var(fid, 'rh',  NF90_DOUBLE, 1, dimids(1), rhvar)
        call wrap_def_var(fid, 'rwet',NF90_DOUBLE, 1, dimids(1), rwetvar)

        dimids(1) = lwdim
        call wrap_def_var(fid, 'lw_band', NF90_DOUBLE, 1, dimids(1), lwvar)

        dimids(1) = swdim
        call wrap_def_var(fid, 'sw_band', NF90_DOUBLE, 1, dimids(1), swvar)

        dimids(1) = wtpdim
        call wrap_def_var(fid, 'wgtpct', NF90_DOUBLE, 1, dimids(1), wtp_var)

        call wrap_put_att_text(fid, rhvar, 'units', 'fraction')
        call wrap_put_att_text(fid, rwetvar, 'units', 'cm')
        call wrap_put_att_text(fid, lwvar, 'units', 'm')
        call wrap_put_att_text(fid, swvar, 'units', 'm')

        call wrap_put_att_text(fid, wtp_var,'units', 'unitless')
        call wrap_put_att_text(fid, wtp_var,'long_name', 'weight percent')

        call wrap_put_att_text(fid, rhvar, 'long_name', 'relative humidity')
        call wrap_put_att_text(fid, rwetvar, 'long_name', 'wet radius')
        call wrap_put_att_text(fid, lwvar, 'long_name', 'longwave bands')
        call wrap_put_att_text(fid, swvar, 'long_name', 'shortwave bands')

        ! Define the variables: abs_lw, ext_sw, ssa_sw, asm_sw
        ! Define 2-dimension (:nrh,:nswbands) LW optics properties: abs_lw, qabs_lw
        dimids(1) = wtpdim
        dimids(2) = lwdim
        call wrap_def_var(fid, 'abs_lw_wtp', NF90_DOUBLE, 2, dimids(1:2), abs_lw_wtp_var)
        call wrap_def_var(fid, 'qabs_lw_wtp',NF90_DOUBLE, 2, dimids(1:2), qabs_lw_wtp_var)

        call wrap_put_att_text(fid, abs_lw_wtp_var, 'units', 'meter^2 kilogram^-1')
        call wrap_put_att_text(fid, qabs_lw_wtp_var,'units', '-')

        ! Define 2-dimension (:nrh,:nswbands) optics properties: ext_sw, qext_sw, ssa_sw, asm_sw
        dimids(1) = wtpdim
        dimids(2) = swdim
        call wrap_def_var(fid, 'ext_sw_wtp', NF90_DOUBLE, 2, dimids(1:2), ext_sw_wtp_var)
        call wrap_def_var(fid, 'qext_sw_wtp',NF90_DOUBLE, 2, dimids(1:2), qext_sw_wtp_var)
        call wrap_def_var(fid, 'ssa_sw_wtp', NF90_DOUBLE, 2, dimids(1:2), ssa_sw_wtp_var)
        call wrap_def_var(fid, 'asm_sw_wtp', NF90_DOUBLE, 2, dimids(1:2), asm_sw_wtp_var)

        call wrap_put_att_text(fid, ssa_sw_wtp_var, 'units', 'fraction')
        call wrap_put_att_text(fid, qext_sw_wtp_var,'units', '-')
        call wrap_put_att_text(fid, ext_sw_wtp_var, 'units', 'meter^2 kilogram^-1')
        call wrap_put_att_text(fid, asm_sw_wtp_var, 'units', '-')

        ! Define the variables for the refractive indicies.
        dimids(1) = swdim
        call wrap_def_var(fid, 'refindex_real_aer_sw', NF90_DOUBLE, 1, dimids(1), sw_r_refidx_var)
        call wrap_def_var(fid, 'refindex_im_aer_sw',   NF90_DOUBLE, 1, dimids(1), sw_i_refidx_var)

        dimids(1) = lwdim
        call wrap_def_var(fid, 'refindex_real_aer_lw', NF90_DOUBLE, 1, dimids(1), lw_r_refidx_var)
        call wrap_def_var(fid, 'refindex_im_aer_lw',   NF90_DOUBLE, 1, dimids(1), lw_i_refidx_var)

        call wrap_put_att_text(fid, sw_r_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, sw_i_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, lw_r_refidx_var, 'units', '-')
        call wrap_put_att_text(fid, lw_i_refidx_var, 'units', '-')

        call wrap_put_att_text(fid, sw_r_refidx_var, 'long_name', 'real refractive index of aerosol - shortwave')
        call wrap_put_att_text(fid, sw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - shortwave')
        call wrap_put_att_text(fid, lw_r_refidx_var, 'long_name', 'real refractive index of aerosol - longwave')
        call wrap_put_att_text(fid, lw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - longwave')


        ! Define fields that define the aerosol properties.
        call wrap_def_dim(fid, 'opticsmethod_len',  32, omdim)
        dimids(1) = omdim
        call wrap_def_var(fid, 'opticsmethod',  NF90_CHAR, 1, dimids(1), omvar)

        call wrap_def_dim(fid, 'namelength',  20, andim)
        dimids(1) = andim
        call wrap_def_var(fid, 'aername',  NF90_CHAR, 1, dimids(1), anvar)

        call wrap_def_dim(fid, 'name_len',  32, namedim)
        dimids(1) = namedim
        call wrap_def_var(fid, 'name',  NF90_CHAR, 1, dimids, namevar)

        call wrap_def_var(fid, 'density',            NF90_DOUBLE, 0, dimids(1), denvar)
        call wrap_def_var(fid, 'sigma_logr',         NF90_DOUBLE, 0, dimids(1), slogvar)
        call wrap_def_var(fid, 'dryrad',             NF90_DOUBLE, 0, dimids(1), dryrvar)
        call wrap_def_var(fid, 'radmin_aer',         NF90_DOUBLE, 0, dimids(1), rminvar)
        call wrap_def_var(fid, 'radmax_aer',         NF90_DOUBLE, 0, dimids(1), rmaxvar)
        call wrap_def_var(fid, 'hygroscopicity',     NF90_DOUBLE, 0, dimids(1), hygrovar)
        call wrap_def_var(fid, 'num_to_mass_ratio',  NF90_DOUBLE, 0, dimids(1), ntmvar)

        call wrap_put_att_text(fid, denvar,   'units', 'kg m^-3')
        call wrap_put_att_text(fid, slogvar,  'units', '-')
        call wrap_put_att_text(fid, dryrvar,  'units', 'm')
        call wrap_put_att_text(fid, rminvar,  'units', 'm')
        call wrap_put_att_text(fid, rmaxvar,  'units', 'm')
        call wrap_put_att_text(fid, hygrovar, 'units', '-')
        call wrap_put_att_text(fid, ntmvar,   'units', 'kg^-1')

        call wrap_put_att_text(fid, denvar,   'long_name', 'aerosol material density')
        call wrap_put_att_text(fid, slogvar,  'long_name', 'geometric standard deviation of aerosol')
        call wrap_put_att_text(fid, dryrvar,  'long_name', 'dry number mode radius of aerosol')
        call wrap_put_att_text(fid, rminvar,  'long_name', 'minimum dry radius of aerosol for bin')
        call wrap_put_att_text(fid, rmaxvar,  'long_name', 'maximum dry radius of aerosol for bin')
        call wrap_put_att_text(fid, hygrovar, 'long_name', 'hygroscopicity of aerosol')
        call wrap_put_att_text(fid, ntmvar,   'long_name', 'ratio of number to mass of aerosol')

        ! End the defintion phase of the netcdf file.
        call wrap_enddef(fid)

        ! Write out the dimensions.
        call wrap_put_var_realx(fid, rhvar, mie_rh(:))
        call wrap_put_var_realx(fid, lwvar, wave(:nlwbands) * 1e-2_f)
        call wrap_put_var_realx(fid, swvar, wave(nlwbands+1:) * 1e-2_f)

        call wrap_put_var_realx(fid, wtp_var, mie_wtp(:)*100._f)

        ! Write out the refractive indicies.
        call wrap_put_var_realx(fid, sw_r_refidx_var, real(refidx(nlwbands+1:, 1)))
        call wrap_put_var_realx(fid, sw_i_refidx_var, aimag(refidx(nlwbands+1:, 1)))
        call wrap_put_var_realx(fid, lw_r_refidx_var, real(refidx(:nlwbands, 1)))
        call wrap_put_var_realx(fid, lw_i_refidx_var, aimag(refidx(:nlwbands, 1)))

        ! Pad the names out with spaces.
        aer_name = '                                '
        aer_name(1:len(trim(c_name))) = c_name

        start_text(1) = 1
        count_text(1) = 32
        call wrap_put_vara_text(fid, namevar, start_text, count_text, (/ aer_name /))
        count_text(1) = 20
        call wrap_put_vara_text(fid, anvar, start_text, count_text, (/ aer_name /))

        count_text(1) = len('hygroscopic_wtp                 ')
        call wrap_put_vara_text(fid, omvar, start_text, count_text, (/ 'hygroscopic_wtp                 ' /))

        call wrap_put_var_realx(fid, denvar,   (/ rho(ibin) * 1e-3_f / 1e-6_f /))
        call wrap_put_var_realx(fid, slogvar,  (/ 0._f /))
        call wrap_put_var_realx(fid, dryrvar,  (/ r(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, rminvar,  (/ rlow(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, rmaxvar,  (/ rup(ibin) * 1e-2_f /))
        call wrap_put_var_realx(fid, hygrovar, (/ 0.6_f /))
        call wrap_put_var_realx(fid, ntmvar,   (/ 1._f / rmass(ibin) / 1e-3_f /))

        ! For now, ext_sw(:nrh, :nswbands) and ext_sw_coreshell(:nrh, :nswbands, :ncoreshellratio) both are calculated
        ! Since other aerosols in CAM may use ext_sw rather than ext_sw_coreshell
        ! Modified by Pengfei Yu
        ! April.1, 2012

        ! calculate qext and ext for pure sulfate dependent on weight percent
        ! ideally qext is based on (wgt,temp,wave), however Beyer et al. (1996) Figure 5
        ! shows sulfate density is roughly 0.006 g/cm3/k, I negelet temp dimension, assuming temp = 270 K
        ! In code, sulfate density is precisely calculated to determine wet raidus
        do iwtp = 1, NMIE_WTP

          ! NOTE: Weight percent is normal a result of the getwetr calculation. To build the
          ! table based upon weight percent, we need to pass in the desired value and a
          ! reference temperature. In that case, the RH is ignored.
          call getwetr(carma, igroup, mie_rh(1), r(ibin), rwet, rho(ibin), rhopwet, rc, wgtpct=mie_wtp(iwtp)*100._f, temp=270._f)
          if (rc < 0) call endrun('carma_CreateOpticsFile::wetr failed.')

          ! Calculate at each wavelength.
          do iwave = 1, NWAVE

            ! Using Mie code, calculate the optical properties: extinction coefficient,
            ! single scattering albedo and asymmetry factor.
            ! Assume the particle is homogeneous (no core).
            !
            ! NOTE: The refractive index for sulfate changes with RH/weight percent, which
            ! is not reflected in this code.
            call mie(carma, &
                     imiertn, &
                     rwet, &
                     wave(iwave), &
                     0._f, &
                     3.0_f, &
                     0.0_f, &
                     1.0_f, &
                     refidx(iwave, 1), &
                     0.0_f, &
                     refidx(iwave, 1), &
                     Qext, &
                     Qsca, &
                     asym, &
                     rc)
            if (rc < 0) call endrun('carma_CreateOpticsFile::mie failed.')

            ! Calculate  the shortwave and longwave properties?
            !
            ! NOTE: miess is in cgs units, but the optics file needs to be in mks
            ! units, so perform the necessary conversions.
            if (iwave <= nlwbands) then

              ! Longwave just needs absorption: abs_lw.
              qabs_lw_wtp(iwtp, iwave) = (Qext - Qsca)                           ! absorption per particle
              abs_lw_wtp (iwtp, iwave) = (Qext - Qsca) * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f)
            else

              ! Shortwave needs extinction, single scattering albedo and asymmetry factor:
              ! ext_sw, ssa_sw and asm_sw.
              qext_sw_wtp(iwtp, iwave - nlwbands) = Qext                             ! extinction per particle
              ext_sw_wtp (iwtp, iwave - nlwbands) = Qext * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f)
              ssa_sw_wtp (iwtp, iwave - nlwbands) = Qsca / Qext
              asm_sw_wtp (iwtp, iwave - nlwbands) = asym
            end if
          end do ! iwave
        end do  ! iwtp

        ! Write out the longwave fields.
        ret = nf90_put_var(fid, abs_lw_wtp_var,  abs_lw_wtp (:, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_SulfateYu: error writing varid =', fid, abs_lw_wtp_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, qabs_lw_wtp_var, qabs_lw_wtp(:, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_SulfateYu: error writing varid =', qabs_lw_wtp_var
           call handle_error(ret)
        end if

        ! Write out the shortwave fields.
        ret = nf90_put_var(fid, ext_sw_wtp_var, ext_sw_wtp (:, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_SulfateYu: error writing varid =', ext_sw_wtp_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, qext_sw_wtp_var,qext_sw_wtp(:, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_SulfateYu: error writing varid =', qext_sw_wtp_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, ssa_sw_wtp_var, ssa_sw_wtp (:, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_SulfateYu: error writing varid =', ssa_sw_wtp_var
           call handle_error(ret)
        end if

        ret = nf90_put_var(fid, asm_sw_wtp_var, asm_sw_wtp (:, :))
        if (ret /= NF90_NOERR) then
           write(iulog,*)'CARMA_CreateOpticsFile_SulfateYu: error writing varid =', asm_sw_wtp_var
           call handle_error(ret)
        end if

        ! Close the file.
        call wrap_close(fid)
      end if
    end do

    return
  end subroutine CARMAMODEL_CreateOpticsFile_SulfateYu


  !! Called at the end of the timestep after all the columns have been processed to
  !! to allow additional diagnostics that have been stored in pbuf to be output.
  !!
  !! NOTE: This is just keeping track of the changes in the interstitial aerosol,
  !! and does not keep track of the aerosol that flows out the top or bottom of the
  !! model or that moves into cloudborne aerosol.
  !!
  !! NOTE: Output occurs a chunk at a time.
  !!
  !!  @version January-2023
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_CalculateCloudborneDiagnostics(carma, state, pbuf, aerclddiag, rc)

    type(carma_type), intent(in)         :: carma        !! the carma object
    type(physics_state), intent(in)      :: state        !! Physics state variables - before pname
    type(physics_buffer_desc), pointer, intent(in)   :: pbuf(:)      !! physics buffer
    real(r8), intent(out)                :: aerclddiag(pcols,MAXCLDAERDIAG) !! the total cloudborne aerosols, supports up to MAXCLDAERDIAG different values
    integer, intent(out)                 :: rc           !! return code, negative indicates failure

    integer                              :: ncols        !! number of columns in the chunk
    integer                              :: icol         !! column index
    integer                              :: ibin         !! bin index
    integer                              :: ienconc      !! concentration element index
    integer                              :: ncore        !! number of cores
    integer                              :: icorelem(NELEM) !! core element index
    real(r8)                             :: mair(pcols,pver)   !! Mass of air column (kg/m2)
    real(r8)                             :: pureso4(pcols,pver) !! Burden pure sulfate (kg/m2)
    real(r8)                             :: mixso4(pcols,pver)  !! Burden mix sulfate (kg/m2)
    real(r8)                             :: bdbc(pcols,pver)    !! Burden BC sulfate (kg/m2)
    real(r8)                             :: bddust(pcols,pver)  !! Burden Dust sulfate (kg/m2)
    real(r8)                             :: bdoc(pcols,pver)    !! Burden OC sulfate (kg/m2)
    real(r8)                             :: bdsalt(pcols,pver)  !! Burden Salt sulfate (kg/m2)
    real(r8)                             :: bdsoa(pcols,pver)   !! Burden SOA sulfate (kg/m2)
    real(r8), pointer, dimension(:,:)    :: mmr                 !! cloudbourne aerosol mmr  (kg/kg)
    character(len=16)                    :: shortname
    character(len=16)                    :: binname
    character(len=16)                    :: concname
    integer                              :: mmr_ndx
    integer                              :: i

    ! Default return code.
    rc = RC_OK

    pureso4(:,:)     = 0._r8
    mixso4(:,:)      = 0._r8
    aerclddiag(:, :) = 0._r8
    bdbc(:, :)       = 0._r8
    bddust(:, :)     = 0._r8
    bdoc(:, :)       = 0._r8
    bdsalt(:, :)     = 0._r8
    bdsoa(:, :)      = 0._r8

    ! Get the air mass in the column
    !
    ! NOTE convert GRAV from cm/s2 to m/s2.
    ncols = state%ncol
    mair(:ncols,:) = state%pdel(:ncols,:) / (GRAV / 100._r8)

    ! For PRSUL, is just the tendency for the concentration element.
    call CARMAGROUP_Get(carma, I_GRP_PRSUL, rc, ienconc=ienconc)
    call CARMAELEMENT_Get(carma, ienconc, rc, shortname=shortname)

    do ibin = 1, nbin

      write(binname, '(A, I2.2)') "CLD"//trim(shortname), ibin
      mmr_ndx = pbuf_get_index(binname)
      call pbuf_get_field(pbuf, mmr_ndx, mmr)

      pureso4(:ncols,:) = pureso4(:ncols,:) + mmr(:ncols,:) * mair(:ncols,:)
    end do

    ! For MXAER, it is the difference in mass between the concentration element
    ! and the sum of the core masses.
    call CARMAGROUP_Get(carma, I_GRP_MXAER, rc, ienconc=ienconc, ncore=ncore, icorelem=icorelem)
    call CARMAELEMENT_Get(carma, ienconc, rc, shortname=concname)

    do ibin = 1, nbin

      write(binname, '(A, I2.2)') "CLD"//trim(concname), ibin
      mmr_ndx = pbuf_get_index(binname)
      call pbuf_get_field(pbuf, mmr_ndx, mmr)

      mixso4(:ncols,:) = mixso4(:ncols,:) + mmr(:ncols,:) * mair(:ncols,:)

      do i = 1, ncore
        call CARMAELEMENT_Get(carma, icorelem(i), rc, shortname=shortname)

        write(binname, '(A, I2.2)') "CLD"//trim(shortname), ibin
        mmr_ndx = pbuf_get_index(binname)
        call pbuf_get_field(pbuf, mmr_ndx, mmr)

        if (shortname .eq. "MXBC") then
          bdbc(:ncols, :) = bdbc(:ncols, :) + mmr(:ncols,:) * mair(:ncols,:)
        else if (shortname .eq. "MXDUST") then
          bddust(:ncols, :) = bddust(:ncols, :) + mmr(:ncols,:) * mair(:ncols,:)
        else if (shortname .eq. "MXOC") then
          bdoc(:ncols, :) = bdoc(:ncols, :) + mmr(:ncols,:) * mair(:ncols,:)
        else if (shortname .eq. "MXSALT") then
          bdsalt(:ncols, :) = bdsalt(:ncols, :) + mmr(:ncols,:) * mair(:ncols,:)
        else if (shortname .eq. "MXSOA") then
          bdsoa(:ncols, :) = bdsoa(:ncols, :) + mmr(:ncols,:) * mair(:ncols,:)
        end if
      end do
    end do

    do icol = 1, ncols
      aerclddiag(icol, 1) = sum(pureso4(icol,:))
      aerclddiag(icol, 2) = sum(mixso4(icol,:))
      aerclddiag(icol, 3) = sum(bdbc(icol,:))
      aerclddiag(icol, 4) = sum(bddust(icol,:))
      aerclddiag(icol, 5) = sum(bdoc(icol,:))
      aerclddiag(icol, 6) = sum(bdsalt(icol,:))
      aerclddiag(icol, 7) = sum(bdsoa(icol,:))
    end do

    return
  end subroutine CARMAMODEL_CalculateCloudborneDiagnostics


  !! Called at the end of the timestep after all the columns have been processed to
  !! to allow additional diagnostics that have been stored in pbuf to be output.
  !!
  !! NOTE: This is just keeping track of the changes in the interstitial aerosol,
  !! and does not keep track of the aerosol that flows out the top or bottom of the
  !! model or that moves into cloudborne aerosol.
  !!
  !! NOTE: Output occurs a chunk at a time.
  !!
  !!  @version January-2023
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_OutputBudgetDiagnostics(carma, icnst4elem, icnst4gas, state, ptend, old_cflux, cflux, dt, pname, rc)
    use cam_history,  only: outfld
    use constituents, only: pcnst, cnst_get_ind

    type(carma_type), intent(in)         :: carma        !! the carma object
    integer, intent(in)                  :: icnst4elem(NELEM, NBIN) !! constituent index for a carma element
    integer, intent(in)                  :: icnst4gas(NGAS)         !! constituent index for a carma gas
    type(physics_state), intent(in)      :: state        !! Physics state variables - before pname
    type(physics_ptend), intent(in)      :: ptend        !! indivdual parameterization tendencies
    real(r8)                             :: old_cflux(pcols,pcnst)  !! cam_in%clfux from before the timestep_tend
    real(r8)                             :: cflux(pcols,pcnst)  !! cam_in%clfux from after the timestep_tend
    real(r8), intent(in)                 :: dt           !! timestep (s)
    character(*), intent(in)             :: pname        !! short name of the physics package
    integer, intent(out)                 :: rc           !! return code, negative indicates failure

    integer                              :: icol         !! column index
    integer                              :: ibin         !! bin index
    integer                              :: i
    integer                              :: icnst        !! constituent index
    integer                              :: ienconc      !! concentration element index
    integer                              :: ncore        !! number of cores
    integer                              :: icorelem(NELEM) !! core element index
    real(r8)                             :: mair(pver)   !! Mass of air column (kg/m2)
    real(r8)                             :: puretend(pcols) !! Tendency pure sulfate (kg/m2/s)
    real(r8)                             :: mixtend(pcols)  !! Tendency mix sulfate (kg/m2/s)
    real(r8)                             :: bdprso4(pcols)  !! Burden pure sulfate (kg/m2)
    real(r8)                             :: bdmxso4(pcols)  !! Burden mixed sulfate (kg/m2)
    real(r8)                             :: cprflux(pcols)  !! Surface Flux tendency, pure sulfate (kg/m2/s)
    real(r8)                             :: cmxflux(pcols)  !! Surface Flux tendency, mix sulfate (kg/m2/s)
    real(r8)                             :: gastend(pcols)  !! Tendency H2SO4 gas (kg/m2/s)
    real(r8)                             :: so2tend(pcols)  !! Tendency SO2 gas (kg/m2/s)
    real(r8)                             :: tottend(pver)   !! Total Tendency mix sulfate (kg/m2/s)

    ! Default return code.
    rc = RC_OK

    puretend(:) = 0._r8
    mixtend(:)  = 0._r8
    gastend(:)  = 0._r8
    so2tend(:)  = 0._r8
    cprflux(:)  = 0._r8
    cmxflux(:)  = 0._r8

    bdmxso4(:)  = 0._r8
    bdprso4(:)  = 0._r8

    ! Add up the sulfate tendencies.
    do icol = 1, state%ncol

      ! Get the air mass in the column
      !
      ! NOTE convert GRAV from cm/s2 to m/s2.
      mair(:) = state%pdel(icol,:) / (GRAV / 100._r8)

      do ibin = 1, nbin

        ! For PRSUL, is just the tendency for the concentration element.
        call CARMAGROUP_Get(carma, I_GRP_PRSUL, rc, ienconc=ienconc)
        icnst = icnst4elem(ienconc, ibin)

        if (ptend%lq(icnst)) then
          puretend(icol) = puretend(icol) + sum(ptend%q(icol,:,icnst) * mair(:))
        end if
        bdprso4(icol) = bdprso4(icol) + sum(state%q(icol,:,icnst) * mair(:))

        cprflux = cprflux(icol) + (cflux(icol,icnst) - old_cflux(icol,icnst))

        ! For MXAER, it is the difference in mass between the concentration element
        ! and the sum of the core masses.
        call CARMAGROUP_Get(carma, I_GRP_MXAER, rc, ienconc=ienconc, ncore=ncore, icorelem=icorelem)
        icnst = icnst4elem(ienconc, ibin)

        tottend(:) = 0._r8
        if (ptend%lq(icnst)) then
          tottend(:) = ptend%q(icol, :, icnst) * mair(:)
        end if
        bdmxso4(icol) = bdmxso4(icol) + sum(state%q(icol,:,icnst) * mair(:))

        cmxflux(icol) = cmxflux(icol) + (cflux(icol,icnst) - old_cflux(icol,icnst))

        do i = 1, ncore
          icnst = icnst4elem(icorelem(i), ibin)
          if (ptend%lq(icnst)) then
            tottend(:) = tottend(:) - ptend%q(icol,:,icnst) * mair(:)
          end if
        end do

        mixtend(icol) = mixtend(icol) + sum(tottend(:))
      end do

      ! Calculate the H2SO4 change.
      icnst = icnst4gas(I_GAS_H2SO4)
      if (ptend%lq(icnst)) then
        gastend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if

      ! Also do SO2
      call cnst_get_ind("SO2", icnst)
      if (ptend%lq(icnst)) then
        so2tend(icol) = sum(ptend%q(icol,:,icnst) * mair(:))
      end if

    end do

    if (carma_do_package_diags) then
       ! Output the total sulfate and H2SO4 tendencies for this physics package.
       call outfld("SO4PRTC_"//trim(pname), puretend(:), pcols, state%lchnk)
       call outfld("SO4MXTC_"//trim(pname), mixtend(:), pcols, state%lchnk)
       call outfld("H2SO4TC_"//trim(pname), gastend(:), pcols, state%lchnk)
       call outfld("SO2TC_"//trim(pname), so2tend(:), pcols, state%lchnk)
       call outfld("SO4PRSF_"//trim(pname), cprflux(:), pcols, state%lchnk)
       call outfld("SO4MXSF_"//trim(pname), cmxflux(:), pcols, state%lchnk)
       call outfld("SO4PRBD_"//trim(pname), bdprso4(:), pcols, state%lchnk)
       call outfld("SO4MXBD_"//trim(pname), bdmxso4(:), pcols, state%lchnk)
    endif

    return
  end subroutine CARMAMODEL_OutputBudgetDiagnostics


  !! Called at the end of the timestep after all the columns have been processed to
  !! to allow additional diagnostics that have been stored in pbuf to be output.
  !!
  !! NOTE: This is just keeping track of the changes in the interstitial aerosol,
  !! and does not keep track of the aerosol that flows out the top or bottom of the
  !! model or that moves into cloudborne aerosol.
  !!
  !! NOTE: Output occurs a chunk at a time.
  !!
  !!  @version January-2023
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_OutputCloudborneDiagnostics(carma, state, pbuf, dt, pname, oldaerclddiag, rc)
    use cam_history, only: outfld

    type(carma_type), intent(in)         :: carma        !! the carma object
    type(physics_state), intent(in)      :: state        !! Physics state variables - before CARMA
    type(physics_buffer_desc), pointer, intent(in)   :: pbuf(:)      !! physics buffer
    real(r8), intent(in)                 :: dt           !! timestep (s)
    character(*), intent(in)             :: pname        !! short name of the physics package
    real(r8), intent(in )                :: oldaerclddiag(pcols,MAXCLDAERDIAG) !! the before timestep cloudborne aerosol diags
    integer, intent(out)                 :: rc           !! return code, negative indicates failure

    real(r8)             :: aerclddiag(pcols,MAXCLDAERDIAG) !! the after timestep cloudborne aerosol diags

    ! Default return code.
    rc = RC_OK

    ! Get the current diagnostics for the cloudborne aerosols.
    call CARMAMODEL_CalculateCloudborneDiagnostics(carma, state, pbuf, aerclddiag, rc)

    ! Output the total sulfate and H2SO4 tendencies for this physics package.
    call outfld("SO4PRCLDTC_"//trim(pname), (aerclddiag(:,1) - oldaerclddiag(:,1)) / dt, pcols, state%lchnk)
    call outfld("SO4MXCLDTC_"//trim(pname), (aerclddiag(:,2) - oldaerclddiag(:,2)) / dt, pcols, state%lchnk)

    ! To be similar to interstitial, where the burden is calculated from the
    ! state before the tendencies are applied, report the old burden not the
    ! current burden.
    ! call outfld("SO4PRCLDBD_"//trim(pname), aerclddiag(:,1), pcols, state%lchnk)
    ! call outfld("SO4MXCLDBD_"//trim(pname), aerclddiag(:,2), pcols, state%lchnk)
    call outfld("SO4PRCLDBD_"//trim(pname), oldaerclddiag(:,1), pcols, state%lchnk)
    call outfld("SO4MXCLDBD_"//trim(pname), oldaerclddiag(:,2), pcols, state%lchnk)

    return
  end subroutine CARMAMODEL_OutputCloudborneDiagnostics


  !! Called at the end of the timestep after all the columns have been processed to
  !! to allow additional diagnostics that have been stored in pbuf to be output.
  !!
  !! NOTE: Output occurs a chunk at a time.
  !!
  !!  @version January-2023
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_OutputDiagnostics(carma, icnst4elem, state, ptend, pbuf, cam_in, rc)
    use cam_history,   only: outfld
    use constituents,  only: cnst_get_ind
    use camsrfexch,    only: cam_in_t

    type(carma_type), intent(in)         :: carma        !! the carma object
    integer, intent(in)                  :: icnst4elem(NELEM, NBIN) !! constituent index for a carma element
    type(physics_state), intent(in)      :: state        !! Physics state variables - before CARMA
    type(physics_ptend), intent(in)      :: ptend        !! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer, intent(in)   :: pbuf(:)  !! physics buffer
    type(cam_in_t), intent(in)           :: cam_in       !! surface inputs
    integer, intent(out)                 :: rc           !! return code, negative indicates failure

    integer                              :: icol         !! column index
    integer                              :: ibin         !! bin index
    real(r8), pointer, dimension(:,:)    :: soacm        !! aerosol tendency due to gas-aerosol exchange  kg/kg/s
    real(r8), pointer, dimension(:,:)    :: soapt        !! aerosol tendency due to no2 photolysis  kg/kg/s
    character(len=16)                    :: binname      !! names bins
    real(r8)                             :: aerclddiag(pcols,MAXCLDAERDIAG) !! the before timestep cloudborne aerosol diags
    integer                              :: i
    integer                              :: icnst        !! constituent index
    integer                              :: ienconc      !! concentration element index
    integer                              :: ncore        !! number of cores
    integer                              :: icorelem(NELEM) !! core element index
    real(r8)                             :: mair(pver)   !! Mass of air column (kg/m2)
    real(r8)                             :: pureso4(pcols) !! pure sulfate (kg/m2)
    real(r8)                             :: mixso4(pcols)  !! mix sulfate (kg/m2)
    real(r8)                             :: cprflux(pcols) !! Surface Flux pure sulfate (kg/m2/s)
    real(r8)                             :: cmxflux(pcols) !! Surface Flux mix sulfate (kg/m2/s)
    real(r8)                             :: h2so4(pcols)   !! H2SO4 gas (kg/m2)
    real(r8)                             :: so2(pcols)     !! SO2 gas (kg/m2)
    real(r8)                             :: bdbc(pcols)    !! Burden BC sulfate (kg/m2)
    real(r8)                             :: bddust(pcols)  !! Burden dust (kg/m2)
    real(r8)                             :: bdoc(pcols)    !! Burden OC sulfate (kg/m2)
    real(r8)                             :: bdsalt(pcols)  !! Burden SALT sulfate (kg/m2)
    real(r8)                             :: bdsoa(pcols)   !! Burden SOA sulfate (kg/m2)
    real(r8)                             :: pureso4mr(pcols,pver) !! Mixing ratio pure sulfate (kg/kg)
    real(r8)                             :: mixso4mr(pcols,pver)  !! Mixing ratio mix sulfate (kg/kg)
    real(r8)                             :: bcmr(pcols,pver)      !! Mixing ratio BC sulfate (kg/kg)
    real(r8)                             :: dustmr(pcols,pver)    !! Mixing ratio dust (kg/kg)
    real(r8)                             :: ocmr(pcols,pver)      !! Mixing ratio OC sulfate (kg/kg)
    real(r8)                             :: saltmr(pcols,pver)    !! Mixing ratio SALT sulfate (kg/kg)
    real(r8)                             :: soamr(pcols,pver)     !! Mixing ratio SOA sulfate (kg/kg)
    character(len=16)                    :: shortname

    ! Default return code.
    rc = RC_OK

    ! Provide diagnostics on the SOA tendencies that affect MXSOA.
    do ibin = 1, NBIN
      write(binname, '(A, I2.2)') "MXSOA", ibin

      call pbuf_get_field(pbuf, ipbuf4soacm(ibin), soacm)
      call outfld(trim(binname)//'CM', soacm(:, :), pcols, state%lchnk)

      call pbuf_get_field(pbuf, ipbuf4soapt(ibin), soapt)
      call outfld(trim(binname)//'PT', soapt(:, :), pcols, state%lchnk)
    end do

    if (carma_do_budget_diags) then
       ! Output the cloudborne SO4 burdens.
       call CARMAMODEL_CalculateCloudborneDiagnostics(carma, state, pbuf, aerclddiag, rc)
       call outfld("SO4PRCLDBD", aerclddiag(:,1), pcols, state%lchnk)
       call outfld("SO4MXCLDBD", aerclddiag(:,2), pcols, state%lchnk)
       call outfld("MXBCCLDBD", aerclddiag(:,3), pcols, state%lchnk)
       call outfld("MXDUSTCLDBD", aerclddiag(:,4), pcols, state%lchnk)
       call outfld("MXOCCLDBD", aerclddiag(:,5), pcols, state%lchnk)
       call outfld("MXSALTCLDBD", aerclddiag(:,6), pcols, state%lchnk)
       call outfld("MXSOACLDBD", aerclddiag(:,7), pcols, state%lchnk)
    endif

    ! Output the interstitial SO4 burdens.
    pureso4(:) = 0._r8
    mixso4(:)  = 0._r8
    cprflux(:) = 0._r8
    cmxflux(:) = 0._r8
    h2so4(:)   = 0._r8
    so2(:)     = 0._r8
    bdbc(:)    = 0._r8
    bddust(:)  = 0._r8
    bdoc(:)    = 0._r8
    bdsalt(:)  = 0._r8
    bdsoa(:)   = 0._r8

    ! Output the mixing ratio
    pureso4mr(:,:) = 0._r8
    mixso4mr(:,:)  = 0._r8
    bcmr(:,:)      = 0._r8
    dustmr(:,:)    = 0._r8
    ocmr(:,:)      = 0._r8
    saltmr(:,:)    = 0._r8
    soamr(:,:)     = 0._r8

    ! Add up the sulfate tendencies.
    do icol = 1, state%ncol

      ! Get the air mass in the column
      !
      ! NOTE convert GRAV from cm/s2 to m/s2.
      mair(:) = state%pdel(icol,:) / (GRAV / 100._r8)

      do ibin = 1, nbin

        ! For PRSUL, is just the tendency for the concentration element.
        call CARMAGROUP_Get(carma, I_GRP_PRSUL, rc, ienconc=ienconc)
        icnst = icnst4elem(ienconc, ibin)

        pureso4mr(icol,:) = pureso4mr(icol,:) + state%q(icol,:,icnst)
        pureso4(icol) = pureso4(icol) + sum(state%q(icol,:,icnst) * mair(:))

        cprflux = cprflux + cam_in%cflx(icol,icnst)

        ! For MXAER, it is the difference in mass between the concentration element
        ! and the sum of the core masses.
        call CARMAGROUP_Get(carma, I_GRP_MXAER, rc, ienconc=ienconc, ncore=ncore, icorelem=icorelem)
        icnst = icnst4elem(ienconc, ibin)

        mixso4mr(icol,:) = mixso4mr(icol,:) + state%q(icol, :, icnst)
        mixso4(icol) = mixso4(icol) + sum(state%q(icol, :, icnst) * mair(:))

        cmxflux(icol) = cmxflux(icol) + cam_in%cflx(icol,icnst)

        do i = 1, ncore
          icnst = icnst4elem(icorelem(i), ibin)

          call CARMAELEMENT_Get(carma, icorelem(i), rc, shortname=shortname)
          if (shortname .eq. "MXBC") then
            bcmr(icol,:) = bcmr(icol,:) + state%q(icol,:,icnst)
            bdbc(icol) = bdbc(icol) + sum(state%q(icol,:,icnst) * mair(:))
          else if (shortname .eq. "MXDUST") then
            dustmr(icol,:) = dustmr(icol,:) + state%q(icol,:,icnst)
            bddust(icol) = bddust(icol) + sum(state%q(icol,:,icnst) * mair(:))
          else if (shortname .eq. "MXOC") then
            ocmr(icol,:) = ocmr(icol,:) + state%q(icol,:,icnst)
            bdoc(icol) = bdoc(icol) + sum(state%q(icol,:,icnst) * mair(:))
          else if (shortname .eq. "MXSALT") then
            saltmr(icol,:) = saltmr(icol,:) + state%q(icol,:,icnst)
            bdsalt(icol) = bdsalt(icol) + sum(state%q(icol,:,icnst) * mair(:))
          else if (shortname .eq. "MXSOA") then
            soamr(icol,:) = soamr(icol,:) + state%q(icol,:,icnst)
            bdsoa(icol) = bdsoa(icol) + sum(state%q(icol,:,icnst) * mair(:))
          end if

        end do
      end do

      ! Calculate the H2SO4 burden.
      call cnst_get_ind("H2SO4", icnst)
      h2so4(icol) = sum(state%q(icol,:,icnst) * mair(:))

      ! Calculate the SO2 burden.
      call cnst_get_ind("SO2", icnst)
      so2(icol) = sum(state%q(icol,:,icnst) * mair(:))
    end do

    if (carma_do_budget_diags) then
       ! Output the total aerosol and gas burdens and the aerosol fluxes.
       call outfld("SO4PRBD", pureso4(:), pcols, state%lchnk)
       call outfld("SO4MXBD", mixso4(:), pcols, state%lchnk)
       call outfld("SO4PRSF", cprflux(:), pcols, state%lchnk)
       call outfld("SO4MXSF", cmxflux(:), pcols, state%lchnk)
       call outfld("H2SO4BD", h2so4(:), pcols, state%lchnk)
       call outfld("SO2BD", so2(:), pcols, state%lchnk)
       call outfld("MXBCBD", bdbc(:), pcols, state%lchnk)
       call outfld("MXDUSTBD", bddust(:), pcols, state%lchnk)
       call outfld("MXOCBD", bdoc(:), pcols, state%lchnk)
       call outfld("MXSALTBD", bdsalt(:), pcols, state%lchnk)
       call outfld("MXSOABD", bdsoa(:), pcols, state%lchnk)
    endif

    ! Output the total aerosol mixing ratio
    call outfld("SO4PRMR", pureso4mr(:,:), pcols, state%lchnk)
    call outfld("MXSO4MR", mixso4mr(:,:), pcols, state%lchnk)
    call outfld("MXBCMR", bcmr(:,:), pcols, state%lchnk)
    call outfld("MXDUSTMR", dustmr(:,:), pcols, state%lchnk)
    call outfld("MXOCMR", ocmr(:,:), pcols, state%lchnk)
    call outfld("MXSALTMR", saltmr(:,:), pcols, state%lchnk)
    call outfld("MXSOAMR", soamr(:,:), pcols, state%lchnk)

    return
  end subroutine CARMAMODEL_OutputDiagnostics



  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011
  !!  @author  Chuck Bardeen
  subroutine CARMAMODEL_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch, only: cam_out_t

    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure

    integer    :: icol

    ! Default return code.
    rc = RC_OK

    return
  end subroutine CARMAMODEL_WetDeposition


  !! Calculates the emissions for CARMA sea salt aerosol particles.
  !!
  !! @author  Tianyi Fan, Chuck Bardeen, Pengfei Yu
  !! @version Dec-2010
  !! originally calculate sea salt flux in EmitParticle, Pengfei Yu make
  !! it a separate subroutine since multiple aerosol types need salt flux
  !! e.g. sea salt, sea salt sulfate, marine organics
  subroutine CARMAMODEL_SaltFlux(carma, ibin, state, r, dr, rmass, cam_in, SaltFlux, rc)
    use ppgrid,        only: pcols
    use physics_types, only: physics_state
    use camsrfexch,    only: cam_in_t

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ibin                  !! bin index
    type(physics_state), intent(in)    :: state                 !! physics state
    real(r8), intent(in)               :: r                     !! bin center (cm)
    real(r8), intent(in)               :: dr                    !! bin width (cm)
    real(r8), intent(in)               :: rmass                 !! bin mass (g)
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: SaltFlux(pcols)       !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index


    ! -------- local variables added for sea salt model ------------
    real(r8)            :: rdrycm, rdry                       ! dry radius [cm], [um]
    real(r8)            :: r80cm, r80                         ! wet radius at relatige humidity of 80% [cm]
    real(r8)            :: ncflx                              ! dF/dr [#/m2/s/um]
    real(r8)            :: Monahan, Clarke, Smith             ! dF/dr [#/m2/s/um]
    real(r8)            :: A_para, B_para, sita_para          ! A, B, and sita parameters in Gong
    real(r8)            :: B_mona                             ! the parameter used in Monahan
    real(r8)            :: W_Caff                             ! Correction factor in Caffrey
    real(r8)            :: u14, ustar_smith, cd_smith         ! 14m wind velocity, friction velocity, and drag coefficient as desired by Andreas source function
    real(r8)            :: wcap                               ! whitecap coverage
    real(r8)            :: fref                               ! correction factor suggested by Hoppe2005
    real(r8), parameter :: xkar = 0.4_r8                      ! Von Karman constant
    real(r8)            :: u10in                              ! 10 meter wind speed use in the emission rate

    ! ------------------------------------------------------------------------------------------------
    ! -- Martensson source function. Coefficients for the parameterization of Ak(c4-c0) and Bk(d4-d0)
    ! -------------------------------------------------------------------------------------------------
    real(r8), parameter :: c41 = -2.576e35_r8
    real(r8), parameter :: c42 = -2.452e33_r8
    real(r8), parameter :: c43 =  1.085e29_r8
    real(r8), parameter :: c31 =  5.932e28_r8
    real(r8), parameter :: c32 =  2.404e27_r8
    real(r8), parameter :: c33 = -9.841e23_r8
    real(r8), parameter :: c21 = -2.867e21_r8
    real(r8), parameter :: c22 = -8.148e20_r8
    real(r8), parameter :: c23 =  3.132e18_r8
    real(r8), parameter :: c11 = -3.003e13_r8
    real(r8), parameter :: c12 =  1.183e14_r8
    real(r8), parameter :: c13 = -4.165e12_r8
    real(r8), parameter :: c01 = -2.881e6_r8
    real(r8), parameter :: c02 = -6.743e6_r8
    real(r8), parameter :: c03 =  2.181e6_r8
    real(r8), parameter :: d41 = 7.188e37_r8
    real(r8), parameter :: d42 = 7.368e35_r8
    real(r8), parameter :: d43 = -2.859e31_r8
    real(r8), parameter :: d31 =-1.616e31_r8
    real(r8), parameter :: d32 =-7.310e29_r8
    real(r8), parameter :: d33 = 2.601e26_r8
    real(r8), parameter :: d21 = 6.791e23_r8
    real(r8), parameter :: d22 = 2.528e23_r8
    real(r8), parameter :: d23 =-8.297e20_r8
    real(r8), parameter :: d11 = 1.829e16_r8
    real(r8), parameter :: d12 =-3.787e16_r8
    real(r8), parameter :: d13 = 1.105e15_r8
    real(r8), parameter :: d01 = 7.609e8_r8
    real(r8), parameter :: d02 = 2.279e9_r8
    real(r8), parameter :: d03 =-5.800e8_r8

    ! ------------------------------------------------------------
    ! ----  Clarke Source Function. Coefficients for Ai    -------
    ! ------------------------------------------------------------
    real(r8), parameter :: beta01 =-5.001e3_r8
    real(r8), parameter :: beta11 = 0.808e6_r8
    real(r8), parameter :: beta21 =-1.980e7_r8
    real(r8), parameter :: beta31 = 2.188e8_r8
    real(r8), parameter :: beta41 =-1.144e9_r8
    real(r8), parameter :: beta51 = 2.290e9_r8
    real(r8), parameter :: beta02 = 3.854e3_r8
    real(r8), parameter :: beta12 = 1.168e4_r8
    real(r8), parameter :: beta22 =-6.572e4_r8
    real(r8), parameter :: beta32 = 1.003e5_r8
    real(r8), parameter :: beta42 =-6.407e4_r8
    real(r8), parameter :: beta52 = 1.493e4_r8
    real(r8), parameter :: beta03 = 4.498e2_r8
    real(r8), parameter :: beta13 = 0.839e3_r8
    real(r8), parameter :: beta23 =-5.394e2_r8
    real(r8), parameter :: beta33 = 1.218e2_r8
    real(r8), parameter :: beta43 =-1.213e1_r8
    real(r8), parameter :: beta53 = 4.514e-1_r8

    ! ---------------------------------------------
    ! coefficient A1, A2 in Andreas's Source funcion
    ! ---------------------------------------------
    real(r8)            ::A1A92
    real(r8)            ::A2A92

    ! ---------------------------------------------
    ! coefficient in Smith's Source funcion
    ! ---------------------------------------------
    real(r8), parameter ::  f1 = 3.1_r8
    real(r8), parameter ::  f2 = 3.3_r8
    real(r8), parameter ::  r1 = 2.1_r8
    real(r8), parameter ::  r2 = 9.2_r8
    real(r8), parameter ::  delta = 10._r8

    ! --------------------------------------------------------------------
    ! ---- constants in calculating the particle wet radius [Gerber, 1985]
    ! --------------------------------------------------------------------
    real(r8), parameter :: c1   = 0.7674_r8        ! .
    real(r8), parameter :: c2   = 3.079_r8         ! .
    real(r8), parameter :: c3   = 2.573e-11_r8     ! .
    real(r8), parameter :: c4   = -1.424_r8        ! constants in calculating the particle wet radius

    ! Default return code.
    rc = RC_OK

    ncol = state%ncol

    ! Add any surface flux here.
    SaltFlux(:ncol) = 0.0_r8

    ! Are we configured for one of the known emission schemes?
    if( carma_seasalt_emis .ne. "Gong"       .and. &
        carma_seasalt_emis .ne. "Martensson" .and. &
        carma_seasalt_emis .ne. "Clarke"     .and. &
        carma_seasalt_emis .ne. "Andreas"    .and. &
        carma_seasalt_emis .ne. "Caffrey"    .and. &
        carma_seasalt_emis .ne. "CMS"        .and. &
        carma_seasalt_emis .ne. "NONE"       .and. &
        carma_seasalt_emis .ne. "CONST"        ) then

       call endrun('carma_EmitParticle:: Invalid sea salt emission scheme.')
    end if

    !**********************************
    ! wet sea salt radius at RH = 80%
    !**********************************
    r80cm   = (c1 *  (r) ** c2 / (c3 * r ** c4 - log10(0.8_r8)) + (r)**3) ** (1._r8/3._r8) ! [cm]
    rdrycm  = r  ! [cm]
    r80     = r80cm *1.e4_r8    ! [um]
    rdry    = rdrycm*1.e4_r8  ! [um]

    do icol = 1,ncol

       ! Only generate sea salt over the ocean.
       if (cam_in%ocnfrac(icol) > 0._r8) then

          !**********************************
          !    WIND for seasalt production
          !**********************************
          call CARMAMODEL_SurfaceWind_salt(icol, cam_in, u10in, rc)

          ! Add any surface flux here.
          ncflx       = 0.0_r8
          Monahan     = 0.0_r8
          Clarke      = 0.0_r8
          Smith       = 0.0_r8

          !**********************************
          !        Whitecap Coverage
          !**********************************
          wcap = 3.84e-6_r8 * u10in ** 3.41_r8      ! in percent, ie., 75%, wcap = 0.75

          !****************************************
          !        Hoppel correction factor
          !        Smith drag coefficients and etc
          !****************************************
          if (u10in .le. 10._r8) then
             cd_smith = 1.14e-3_r8
          else
             cd_smith = (0.49_r8 + 0.065_r8 * u10in) * 1.e-3_r8
          end if

          ! ustar_smith = cd_smith **0.5_r8 * u10in
          !
          ! We don't have vg yet, since that is calculated by CARMA. That will require
          ! a different interface for the emissions, storing vg in the physics buffer,
          ! and/or doing some duplicate calculations for vg assuming 80% RH.
          !          fref = (delta/state%zm(icol, pver))**(vg(icol, ibin, igelem(i))/(xkar*ustar_smith))
          fref = 1.0_r8

          !**********************************
          !        Source Functions
          !**********************************
          if (carma_seasalt_emis .eq. 'NONE') then
             ncflx = 0._r8
          end if

          if (carma_seasalt_emis .eq. 'CONST') then
             ncflx = 1.e-5_r8
          end if

          !-------Gong source function------
          if (carma_seasalt_emis == "Gong") then
             sita_para = 30
             A_para = - 4.7_r8 * (1+ sita_para * r80) ** (- 0.017_r8 * r80** (-1.44_r8))
             B_para = (0.433_r8 - log10(r80)) / 0.433_r8
             ncflx = 1.373_r8* u10in ** 3.41_r8 * r80 ** A_para * (1._r8 + 0.057_r8 * r80**3.45_r8) * 10._r8 ** (1.607_r8 * exp(- B_para **2))
             !            if (do_print) write(LUNOPRT, *) "Gong: ncflx = ", ncflx, ", u10n = ", u10in
          end if

          !------Martensson source function-----
          if (carma_seasalt_emis == "Martensson") then
             if (rdry .le. 0.0725_r8) then
                ncflx = (Ak1(rdry*1.0e-6_r8)* (25._r8+273._r8) + Bk1(rdry*1.0e-6_r8)) * wcap      ! dF/dlogr [#/s/m2]
                ncflx = ncflx / (2.30258509_r8 * rdry)                                            ! dF/dr    [#/s/m2/um]
             elseif (rdry .gt. 0.0725_r8 .and. rdry .le. 0.2095_r8) then
                ncflx = (Ak2(rdry*1.0e-6_r8)* (25._r8+273._r8) + Bk2(rdry*1.0e-6_r8)) * wcap      ! dF/dlogr [#/s/m2]
                ncflx = ncflx / (2.30258509_r8 * rdry)                                            ! dF/dr    [#/s/m2/um]
             elseif (rdry .gt. 0.2095_r8 .and. rdry .le. 1.4_r8) then
                ncflx = (Ak3(rdry*1.0e-6_r8)* (25._r8+273._r8) + Bk3(rdry*1.0e-6_r8)) * wcap      ! dF/dlogr [#/s/m2]
                ncflx = ncflx / (2.30258509_r8 * rdry)                                            ! dF/dr    [#/s/m2/um]
             else
                ncflx = 0._r8
             end if
          end if

          !-------Clarke source function-------
          if (carma_seasalt_emis == "Clarke")then
             if (rdry .lt. 0.066_r8) then
                ncflx = A1(rdry) * 1.e4_r8 * wcap                              ! dF/dlogr [#/s/m2]
                ncflx = ncflx / (2.30258509_r8 * rdry)                        ! dF/dr    [#/s/m2/um]
             elseif (rdry .ge. 0.066_r8 .and. rdry .lt. 0.6_r8) then
                ncflx = A2(rdry) * 1.e4_r8 * wcap                             ! dF/dlogr [#/s/m2]
                ncflx = ncflx / (2.30258509_r8 * rdry)                        ! dF/dr    [#/s/m2/um]
             elseif (rdry .ge. 0.6_r8 .and. rdry .lt. 4.0_r8) then
                ncflx = A3(rdry) * 1.e4_r8 * wcap                             ! dF/dlogr [#/s/m2]
                ncflx= ncflx / (2.30258509_r8 * rdry)                         ! dF/dr    [#/s/m2/um]
             else
                ncflx = 0._r8
             end if
          end if

          !-----------Caffrey source function------------
          if (carma_seasalt_emis == "Caffrey") then

             !Monahan
             B_mona = (0.38_r8 - log10(r80)) / 0.65_r8
             Monahan = 1.373_r8 * (u10in**3.41_r8) * r80**(-3._r8) * (1._r8 + 0.057_r8 *r80**1.05_r8)  * 10._r8 ** (1.19_r8 * exp(-1._r8 * B_mona**2)) ! dF/dr

             !Smith
             u14 = u10in * (1._r8 + cd_smith**0.5_r8 / xkar * log(14._r8 / 10._r8))  ! 14 meter wind
             A1A92 = 10._r8 ** (0.0676_r8 * u14 + 2.430_r8)
             A2A92 = 10._r8 ** (0.9590_r8 * u14**0.5_r8 - 1.476_r8)
             Smith = A1A92*exp(-f1 *(log(r80/r1))**2) + A2A92*exp(-f2 * (log(r80/r2))**2)     ! dF/dr   [#/m2/s/um]

             !Caffrey based on Monahan and Smith
             W_Caff = 1.136_r8 **(-1._r8 * rdry ** (-0.855_r8))*(1._r8 + 0.2_r8/rdry)
             if (rdry .lt. 0.15_r8) then
                ncflx = Monahan
             else
                if (u10in .le. 9._r8) then
                   ncflx = Monahan
                else
                   if(Monahan .ge. Smith) then
                      ncflx = Monahan
                   else
                      ncflx = Smith
                   end if
                end if
             end if

             ncflx = ncflx * W_Caff

             !%%%%%%%%%%%%%%%%%%%%%%%%%
             ! Apply Hoppel correction
             !%%%%%%%%%%%%%%%%%%%%%%%%%
             ncflx = ncflx * fref
          end if

          !--------CMS (Clarke, Monahan, and Smith source function)-------
          if (carma_seasalt_emis == "CMS") then

             !Clarke
             if (rdry .lt. 0.066_r8) then
                Clarke = A1(rdry) * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2]
                Clarke = Clarke / (2.30258509_r8 * rdry)               ! dF/dr    [#/s/m2/um]
             elseif ((rdry .ge. 0.066_r8) .and. (rdry .lt. 0.6_r8)) then
                Clarke = A2(rdry) * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2]
                Clarke = Clarke / (2.30258509_r8 * rdry)               ! dF/dr    [#/s/m2/um]
             elseif ((rdry .ge. 0.6_r8) .and. (rdry .lt. 4.0_r8)) then
                Clarke = A3(rdry) * 1.e4_r8 * wcap                      ! dF/dlogr [#/s/m2]
                Clarke= Clarke / (2.30258509_r8 * rdry)                 ! dF/dr    [#/s/m2/um]
             end if

             !Monahan
             B_Mona = (0.38_r8 - log10(r80)) / 0.65_r8
             Monahan = 1.373_r8 * u10in ** 3.41_r8 * r80 ** (-3._r8) * (1._r8 + 0.057_r8 * r80**1.05_r8) * 10._r8 ** (1.19_r8 * exp(- B_Mona **2))

             !Smith
             u14 = u10in * (1._r8 + cd_smith**0.5_r8 / xkar*log(14._r8 / 10._r8))  ! 14 meter wind
             A1A92 = 10._r8 ** (0.0676_r8 * u14 + 2.430_r8)
             A2A92 = 10._r8 ** (0.9590_r8 * u14**0.5_r8 - 1.476_r8)
             Smith = A1A92*exp(-f1 *(log(r80 / r1))**2) + A2A92*exp(-f2 * (log(r80 / r2))**2)     ! dF/dr   [#/m2/s/um]

             !%%%%%%%%%%%%%%%%%%%%%%%%%
             !     CMS1 or CMS2
             !%%%%%%%%%%%%%%%%%%%%%%%%%
             !          if (rdry .lt. 0.1_r8) then   ! originally cut at 0.1 um
             ! ***CMS1*****
             if (rdry .lt. 1._r8) then    ! cut at 1.0 um
                ! ***CMS2*****
                !          if (rdry .lt. 2._r8) then    ! cut at 2.0 um
                ncflx = Clarke
             else
                if (u10in .lt. 9._r8) then
                   ncflx = Monahan
                else
                   if (Monahan .gt. Smith) then
                      ncflx = Monahan
                   else
                      ncflx = Smith
                   end if
                end if
             end if

             !%%%%%%%%%%%%%%%%%%%%%%%%%
             ! Apply Hoppel correction
             !%%%%%%%%%%%%%%%%%%%%%%%%%
             ncflx = ncflx * fref
          end if

          ! convert ncflx [#/m^2/s/um] to surfaceFlx [kg/m^2/s]
          SaltFlux(icol) = ncflx * dr * rmass * 10._r8      ! *1e4[um/cm] * 1.e-3[kg/g]

          !          if (do_print) write(LUNOPRT, *) "ibin = ", ibin, ", igroup = ", igroup
          !          if (do_print) write(LUNOPRT, *) "dr = ", dr, ", rmass = ", rmass
          !          if (do_print) write(LUNOPRT, *) "ncflx = " , ncflx, ", SaltFlux = ", SaltFlux(icol)

          ! weighted by the ocean fraction
          SaltFlux(icol) = SaltFlux(icol) * cam_in%ocnfrac(icol)
       end if
    end do

  contains

    ! Coefficient Ak in Martensson's source functions
    pure real(r8) function Ak1(rpdry)
      real(r8),intent(in) :: rpdry
      Ak1 = c41*(2._r8*rpdry)**4 + c31*(2._r8*rpdry) ** 3 + c21*(2._r8*rpdry)**2 + c11*(2._r8*rpdry)+ c01
    end function Ak1

    pure real(r8) function Ak2(rpdry)
      real(r8),intent(in) :: rpdry
      Ak2 = c42*(2._r8*rpdry)**4 + c32*(2._r8*rpdry) ** 3 + c22*(2._r8*rpdry)**2 + c12*(2._r8*rpdry)+ c02
    end function Ak2

    pure real(r8) function Ak3(rpdry)
      real(r8),intent(in) :: rpdry
      Ak3 = c43*(2._r8*rpdry)**4 + c33*(2._r8*rpdry) ** 3 + c23*(2._r8*rpdry)**2 + c13*(2._r8*rpdry)+ c03
    end function Ak3

    ! Coefficient Bk in Martensson's source functions
    pure real(r8) function Bk1(rpdry)
      real(r8),intent(in) :: rpdry
      Bk1= d41*(2._r8*rpdry)**4 + d31*(2._r8*rpdry) ** 3 + d21*(2._r8*rpdry)**2 + d11*(2._r8*rpdry)+ d01
    end function Bk1

    pure real(r8) function Bk2(rpdry)
      real(r8),intent(in) :: rpdry
      Bk2 = d42*(2._r8*rpdry)**4 + d32*(2._r8*rpdry) ** 3 + d22*(2._r8*rpdry)**2 + d12*(2._r8*rpdry)+ d02
    end function Bk2

    pure real(r8) function Bk3(rpdry)
      real(r8),intent(in) :: rpdry
      Bk3 = d43*(2._r8*rpdry)**4 + d33*(2._r8*rpdry) ** 3 + d23*(2._r8*rpdry)**2 + d13*(2._r8*rpdry)+ d03
    end function Bk3

    ! Coefficient Ak in Clarkes's source function
    pure real(r8) function A1(rpdry)
      real(r8),intent(in) :: rpdry
      A1 = beta01 + beta11*(2._r8*rpdry) + beta21*(2._r8*rpdry)**2 + beta31*(2._r8*rpdry)**3 &
           + beta41*(2._r8*rpdry)**4 + beta51*(2._r8*rpdry)**5
    end function A1

    pure real(r8) function A2(rpdry)
      real(r8),intent(in) :: rpdry
      A2 = beta02 + beta12*(2._r8*rpdry) + beta22*(2._r8*rpdry)**2 + beta32*(2._r8*rpdry)**3 &
           + beta42*(2._r8*rpdry)**4 + beta52*(2._r8*rpdry)**5
    end function A2

    pure real(r8) function A3(rpdry)
      real(r8),intent(in) :: rpdry
      A3 = beta03 + beta13*(2._r8*rpdry) + beta23*(2._r8*rpdry)**2 + beta33*(2._r8*rpdry)**3 &
           + beta43*(2._r8*rpdry)**4 + beta53*(2._r8*rpdry)**5
    end function A3

  end subroutine CARMAMODEL_SaltFlux


  !! Calculate the sea surface wind with a Weibull distribution.
  !!
  !! @author  Tianyi Fan
  !! @version August-2010
  subroutine CARMAMODEL_SurfaceWind_salt(icol, cam_in, u10in, rc)
    use camsrfexch, only: cam_in_t

    ! in and out field
    integer, intent(in)                 :: icol                  !! column index
    type(cam_in_t), intent(in)          :: cam_in                !! surface inputs
    real(r8), intent(out)               :: u10in                 !! the 10m wind speed put into the source function
    integer, intent(out)                :: rc                    !! return code, negative indicates failure

    ! local variables
    real(r8) :: uWB341              ! the nth mean wind with integration using Weibull Distribution(integrate from threshold wind velocity)

    rc = RC_OK

    uWB341 = 0._r8

    ! calc. the Weibull wind distribution
    u10in = cam_in%u10(icol)

    call CARMAMODEL_WeibullWind(u10in, uth_salt, 3.41_r8, uWB341)

    u10in = uWB341 ** (1._r8 / 3.41_r8)

!    if (do_print) write(LUNOPRT, *) 'CARMA_SurfaceWind: icol ',icol, ', u10 =', cam_in%u10(icol), ', u10in =', u10in

    return
  end subroutine CARMAMODEL_SurfaceWind_salt



  !!  Determines the mass fraction for the clay (submicron) bins based upon
  !!  Tegen and Lacis [1996]. The total fraction for all clay bins should
  !!  add up to 1.
  !!
  !!  NOTE: WOuld it be better to interpolate this into the bins rather than
  !!  assigning all CARMA bins within a Tegen & Lacis bin the same value?
  !!
  !!  NOTE: Should any mass go to bins smaller than the smallest one used by
  !!  Tegen and Lacis?
  !!
  !!  @version July-2012
  !!  @author  Lin Su, Pengfei Yu, Chuck Bardeen
  subroutine CARMAMODEL_ClayMassFraction(carma, igroup, rdust, rc)

    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: igroup      !! the carma group index
    real(r8), intent(in)                 :: rdust(NBIN) !! radius assuming entire particle is dust
    integer, intent(inout)               :: rc          !! return code, negative indicates failure

    ! Bins and mass fraction from Tegen and Lacis.
    integer, parameter  :: NBIN_TEGEN = 4
    real(r8)            :: tl_rmin(NBIN_TEGEN) = (/ 1.e-5_r8,  1.8e-5_r8, 3.e-5_r8, 6.e-5_r8 /)
    real(r8)            :: tl_rmax(NBIN_TEGEN) = (/ 1.8e-5_r8, 3.e-5_r8,  6.e-5_r8, 1.e-4_r8 /)
    real(r8)            :: tl_mf(NBIN_TEGEN)   = (/ 0.009_r8,  0.081_r8,  0.234_r8, 0.676_r8 /)

    ! Local Variables
    integer, parameter  :: IBELOW = 1
    integer, parameter  :: IABOVE = 6
    integer             :: tl_count(NBIN_TEGEN+2)  ! count number in Tegen and Lacis ranges
    integer             :: ind_up(NBIN_TEGEN+2)
    integer             :: ind_low(NBIN_TEGEN+2)
    integer             :: j                    ! local index number
    integer             :: ibin                 ! carma bin index

    ! Default return code.
    rc = RC_OK

    ! Figure out how many of the CARMA bins are in each of the Tegen and Lacis
    ! ranges.
    tl_count(:) = 0

    do ibin = 1, NBIN

      ! Smaller than the range.
      if (rdust(ibin) < tl_rmin(1)) then
        tl_count(IBELOW) = tl_count(IBELOW) + 1
      end if

      ! In the range
      do j = 1, NBIN_TEGEN
        if (rdust(ibin) < tl_rmax(j) .and. rdust(ibin) >= tl_rmin(j)) then
          tl_count(j+1) = tl_count(j+1) + 1
        end if
      end do

      ! Bigger than the range.
      if (rdust(ibin) >= tl_rmax(NBIN_TEGEN)) then
        tl_count(IABOVE) = tl_count(IABOVE) + 1
      end if
    end do

    ! Determine where the boundaries are between the TEGEN bins and
    ! the CARMA bin structure.
    ind_up(:)   = 0
    ind_low(:)  = 0
    ind_up (IBELOW)  = tl_count(IBELOW)
    ind_low(IBELOW)  = min(1, tl_count(IBELOW))

    do j = 1, 5
      ind_up (j+1) = ind_up(j) + tl_count(j+1)
      ind_low(j+1) = ind_up(j) + min(tl_count(j+1), 1)
    end do

    ! No mass to bins smaller than the smallest size.
    clay_mf(:) = 0._r8

    ! NOTE: This won't work right if the dust bins are coarser than
    ! the Tegen and Lacis bins. In this case mass fraction would need
    ! to be combined from the Tegen & Lacis bins into a CARMA bin.
    do j = 1, NBIN_TEGEN
      if (tl_count(j+1) > 0) then
        clay_mf(ind_low(j+1):ind_up(j+1)) = tl_mf(j) / tl_count(j+1)
      end if
    end do

    clay_mf(ind_low(IABOVE):) = 1._r8

    return
  end subroutine CARMAMODEL_ClayMassFraction


  !! Calculate the sea surface wind with a Weibull distribution.
  !!
  !! NOTE: This should be combined with a similar routine in the sea salt
  !! model, and any differences should be control by parameters into this
  !! routine (and perhaps namelist variables).
  !!
  !! @author  Lin Su, Pengfei Yu, Chuck Bardeen
  !! @version July-2012
  subroutine CARMAMODEL_SurfaceWind(carma, icol, ielem, igroup, ibin, cam_in, uv10, wwd, uth, rc)
    use camsrfexch, only: cam_in_t

    ! in and out field
    type(carma_type), intent(in)        :: carma                 !! the carma object
    integer, intent(in)                 :: icol                  !! column index
    integer, intent(in)                 :: ielem                 !! element index
    integer, intent(in)                 :: igroup                !! group index
    integer, intent(in)                 :: ibin                  !! bin index
    type(cam_in_t), intent(in)          :: cam_in                !! surface inputs
    real(r8), intent(out)               :: uv10                  !! the 10m wind speed (m/s)
    real(r8), intent(out)               :: wwd                   !! the 10m wind speed  with Weibull applied (m/s)
    real(r8), intent(out)               :: uth                   !! the 10m wind threshold (m/s)
    integer,  intent(inout)             :: rc                    !! return code, negative indicates failure

    real(r8), parameter                 :: vk = 0.4_r8           ! von Karman constant
    real(r8)                            :: rmass(NBIN)           ! CARMA bin mass (g)
    real(r8)                            :: r                     ! CARMA bin center (cm)
    real(r8)                            :: rhop(NBIN)            ! CARMA partile element density (g/cm3)
    real(r8)                            :: uthfact               !
    real(r8), parameter                 :: rhoa = 1.25e-3_r8     ! Air density at surface

    rc = RC_OK

    ! Get the 10 meter wind speed
    uv10 = cam_in%u10(icol)

    ! Calculate the threshold wind speed of each bin [Marticorena and Bergametti,1995]
    ! note that in cgs units --> m/s
    call CARMAGROUP_GET(carma, igroup, rc, rmass=rmass)
    if (RC < RC_ERROR) return

    ! Define particle # concentration element index for current group
    call CARMAELEMENT_Get(carma, ielem, rc, rho=rhop)
    if (RC < RC_ERROR) return

    ! Calculate the radius assuming that all the mass will be emitted as this
    ! element.
    r = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin))**(1._r8 / 3._r8)

    if (cam_in%soilw(icol) >= 0._r8 .AND. cam_in%soilw(icol) < 0.5_r8) then

       ! Prevent small values of soilw from driving uthfact negative, but allow
       ! for dust emissions even when soilw is 0.
       uthfact = 1.2_r8 + 0.2_r8*log10(max(0.001_r8, cam_in%soilw(icol)))

       if (r > 2.825e-5_r8) then  ! r(4) = 2.825e-5 cm
           uth = uthfact * 1.e-2_r8 * 0.13_r8 * sqrt(rhop(ibin)*GRAV*r*2._r8/rhoa) &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/(r*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*(r*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       else
           uth = uthfact*1.e-2_r8* 0.13_r8 * sqrt(rhop(ibin)*GRAV*(.75e-4_r8)*2._r8/rhoa)   &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/((.75e-4_r8)*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*((.75e-4_r8)*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       endif
    else
       uth = uv10
    endif

    ! Use Weibull with Lansing's estimate for shape.
    call CARMAMODEL_WeibullWind(uv10, uth, 2._r8, wwd)

    ! Set the threshold to the weibull wind value if sol moisture >= 0.5,
    ! to turn off emissions.
    if (cam_in%soilw(icol) >= 0.5_r8) then
      uth = sqrt(wwd)
    end if

    return
  end subroutine CARMAMODEL_SurfaceWind


  !! Read in the dust source (soil) erodibility factor from a NETCDF file. In this
  !! processes, the data is regridded from the source size to the size needed by the
  !! model.
  !!
  !! NOTE: This is currently doing 2-D interpolation, but it really should be doing
  !! regridding.
  !!
  !! @author  Pengfei Yu
  !! @version July-2012

!! st
!! could use /components/cam/src/chemistry/aerosol/soil_erod_mod.F90 here insted of this routine?
  subroutine CARMAMODEL_ReadSoilErosionFactor(rc)
    use ppgrid,             only: begchunk, endchunk, pcols
    use ioFileMod,          only: getfil
    use interpolate_data,   only: lininterp_init, lininterp, interp_type, lininterp_finish
    use phys_grid,          only: get_rlon_all_p, get_rlat_all_p, get_ncols_p
    use wrap_nf

    integer, intent(out)                      :: rc                    !! return code, negative indicates failure

    ! local variables
    integer                                   :: idvar, f_nlon, f_nlat, idlat, idlon
    integer                                   :: fid, fid_lon, fid_lat
    real(r8), allocatable, dimension(:,:)     :: ero_factor
    character(len=256)                        :: ero_file
    real(r8), allocatable, dimension(:)       :: ero_lat               ! latitude dimension
    real(r8), allocatable, dimension(:)       :: ero_lon               ! latitude dimension
    type (interp_type)                        :: lat_wght, lon_wght
    real(r8)                                  :: lat(pcols)            ! latitude index
    real(r8)                                  :: lon(pcols)            ! longitude index
    integer                                   :: i
    integer                                   :: lchnk                 ! chunk identifier
    integer                                   :: ncol                  ! number of columns in chunk

    real(r8), parameter   :: zero=0_r8, twopi=2_r8*pi, degs2rads = pi/180._r8

    rc = RC_OK

    ! Open the netcdf file (read only)
    call getfil(carma_soilerosion_file, ero_file, 0)
    call wrap_open(ero_file, 0, fid)

    ! Get file dimensions
    call wrap_inq_dimid(fid, 'plon', fid_lon)
    call wrap_inq_dimid(fid, 'plat', fid_lat)
    call wrap_inq_dimlen(fid, fid_lon, f_nlon)
    call wrap_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(ero_lat(f_nlat))
    allocate(ero_lon(f_nlon))
    allocate(ero_factor (f_nlon, f_nlat))
    allocate(soil_factor(pcols, begchunk:endchunk))

    ! Read in the tables.
    call wrap_inq_varid(fid, 'new_source', idvar)
    i = nf90_get_var (fid, idvar, ero_factor)
    if (i/=NF90_NOERR) then
       write(iulog,*)'CARMA_ReadSoilErosionFactor: error reading varid =', idvar
       call handle_error (i)
    end if
    call wrap_inq_varid(fid, 'plat', idlat)
    call wrap_get_var_realx(fid, idlat,  ero_lat)
    call wrap_inq_varid(fid, 'plon', idlon)
    call wrap_get_var_realx(fid, idlon,  ero_lon)

    ero_lat(:) = ero_lat(:)*degs2rads
    ero_lon(:) = ero_lon(:)*degs2rads

    ! Close the file.
    call wrap_close(fid)

    do lchnk=begchunk, endchunk
       ncol = get_ncols_p(lchnk)

       call get_rlat_all_p(lchnk, pcols, lat)
       call get_rlon_all_p(lchnk, pcols, lon)

       call lininterp_init(ero_lon, f_nlon, lon, ncol, 2, lon_wght, zero, twopi)
       call lininterp_init(ero_lat, f_nlat, lat, ncol, 1, lat_wght)

       call lininterp(ero_factor, f_nlon, f_nlat, soil_factor(1:ncol,lchnk), ncol, lon_wght, lat_wght)

       call lininterp_finish(lon_wght)
       call lininterp_finish(lat_wght)
    end do

    deallocate(ero_lat)
    deallocate(ero_lon)
    deallocate(ero_factor)

  end subroutine CARMAMODEL_ReadSoilErosionFactor

  !! Calculate the nth mean of u using Weibull wind distribution
  !! considering the threshold wind velocity. This algorithm
  !! integrates from uth to infinite (u^n P(u)du )
  !!
  !! @author  Tianyi Fan
  !! @version August-2010
   subroutine CARMAMODEL_WeibullWind(u, uth, n, uwb, wbk)
    use shr_spfn_mod, only: gamma =>  shr_spfn_gamma, igamma => shr_spfn_igamma

    real(r8), intent(in)  :: u      ! mean wind speed
    real(r8), intent(in)  :: uth    ! threshold velocity
    real(r8), intent(in)  :: n      ! the rank of u in the integration
    real(r8), intent(out) :: uwb    ! the Weibull distribution
    real(r8), intent(in), optional ::  wbk    ! the shape parameter

    ! local variable
    real(r8)  :: k                  ! the shape parameter in Weibull distribution
    real(r8)  :: c                  ! the scale parameter in Weibull distribution

    if (present(wbk)) then
      k = wbk
    else
      k = 0.94_r8*u**0.5_r8        ! follow Grini and Zender, 2004JGR
 !    k = 2.5_r8                   ! Lansing's estimate
    end if

    ! If u is 0, then k can be 0, which makes a lot of this undefined.
    ! Just return 0. in this case.
    if (u < 0.35_r8) then
      uwb = 0._r8
    else
      c   = u * (gamma(1._r8 + 1._r8 / k))**(-1._r8)
      uwb = c**n * igamma(n / k + 1._r8, (uth / c)**k)
    end if

  end subroutine CARMAMODEL_WeibullWind

  !! Read BC data from three components:
  !! 1. GAINS anthropogenic; 2. Ship Emission; 3. GFEDv3; 4. Aircraft
  !! GAINS unit: kt/year; 2D; lon:-180-180
  !! Ship Emission unit: kg/m2/s; 3D (month,lat,lon); lon:0-360
  !! GFEDv3 unit: g/m2/month; 3D (month,lat,lon); lon:-180-180
  !!
  !! @author  Pengfei Yu
  !! @version May-2013
  subroutine CARMAMODEL_BCOCRead(rc)
    use pmgrid,        only: plat, plon
    use ioFileMod,     only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
    use pio,            only : file_desc_t, var_desc_t, &
                               pio_inq_dimid, pio_inq_varid, &
                               pio_get_var, pio_nowrite, pio_inq_dimlen, &
                               pio_inq_dimlen, pio_closefile
    use dycore,        only: dycore_is

    integer, intent(out)                      :: rc                    !! return code, negative indicates failure

    ! local variables
    integer                                   :: f_nlon, f_nlat, f_ntime
    integer                                   :: fid_lon, fid_lat, fid_time
    real(r8), allocatable, dimension(:,:)     :: BC_f2d, BC2d, OC_f2d, OC2d
    real(r8), allocatable, dimension(:,:,:)   :: BC_f3d, BC3d, OC_f3d, OC3d
!
    character(len=256)                        :: BC_GAINS_file
    character(len=256)                        :: OC_GAINS_file
    character(len=256)                        :: BC_GFEDv3_file
    character(len=256)                        :: OC_GFEDv3_file
    character(len=256)                        :: BC_ship_file
    character(len=256)                        :: OC_ship_file
!
    real(r8), allocatable, dimension(:,:,:)       :: BC_anthro_GAINS
    real(r8), allocatable, dimension(:,:,:)       :: OC_anthro_GAINS
    real(r8), allocatable, dimension(:,:,:)       :: BC_GFEDv3
    real(r8), allocatable, dimension(:,:,:)       :: OC_GFEDv3
    real(r8), allocatable, dimension(:,:,:)       :: BC_ship_GAINS
    real(r8), allocatable, dimension(:,:,:)       :: OC_ship_GAINS
!
    real(r8), allocatable, dimension(:)       :: BC_lat, OC_lat       ! latitude dimension
    real(r8), allocatable, dimension(:)       :: BC_lon, OC_lon       ! latitude dimension
    type (interp_type)                        :: wgt1, wgt2
    real(r8)                                  :: lat(plat), lon(plon)
    integer                                   :: i, itime
    real(r8)                                                              :: rearth, gridarea
    integer                                                                       :: nmonth
    real(r8)                                                              :: tempor(plon,plat)
    real(r8), allocatable, dimension(:,:,:)       :: tempor3d
    real(r8), allocatable, dimension(:,:)         :: tempor2d
    real(r8), allocatable, dimension(:)           :: tempor1d
    integer                                                                       :: mid_idx
    real(r8), allocatable, dimension(:,:)         :: BC_dom_f2d, OC_dom_f2d
    real(r8), allocatable, dimension(:,:,:)       :: BC_dom_f3d, OC_dom_f3d
    real(r8), allocatable, dimension(:,:,:)       :: BC_awb_f3d, OC_awb_f3d
    real(r8), allocatable, dimension(:,:)         :: BC2d_dom, OC2d_dom
    real(r8), allocatable, dimension(:)           :: facH, facL
    integer                                                                       :: ind_15N, ind_45N, ierr
    type(file_desc_t) :: fid
    type(var_desc_t) :: idvar, idlat, idlon, idvar_dom, idvar_awb

    real(r8) :: nlats

    rc = RC_OK

    if(dycore_is('UNSTRUCTURED') ) then
       call endrun('CARMAMODEL_BCOCRead: Yu2015 emissions not implemented for unstructured grids' )
    end if

    ! get model lat and lon
    nlats = plat-1 ! gnu compiler workaround
    do i = 1, plat
       lat(i) = 180._r8/(nlats)*(i-1)-90._r8
    end do
    do i = 1, plon
       lon(i) = 360._r8/plon*(i-1)
    end do

!
    nmonth = 12

    if(carma_BCOCemissions == 'Yu2015')then
       ! allocate BCnew and OCnew, unit is #/cm2/s
       allocate(BCnew(plat, plon, nmonth))
       allocate(OCnew(plat, plon, nmonth))
       BCnew = -huge(1._r8)
       OCnew = -huge(1._r8)
    endif

! monthly fraction of domestic emission
    allocate(facH(nmonth))
    allocate(facL(nmonth))
    facH = (/0.18_r8,0.14_r8,0.13_r8,0.08_r8,0.04_r8,0.02_r8,0.01_r8,&
                0.02_r8,0.03_r8,0.07_r8,0.11_r8,0.17_r8/)
    facL = (/0.17_r8,0.14_r8,0.11_r8,0.06_r8,0.04_r8,0.04_r8,0.04_r8,&
                0.04_r8,0.04_r8,0.06_r8,0.10_r8,0.15_r8/)

! find index for 15N and 45N
    do i = 1, plat
       if (lat(i) .gt. 15._r8) then
          ind_15N = i
          exit
       endif
    end do
!
    do i = 1, plat
       if (lat(i) .gt. 45._r8) then
          ind_45N = i
          exit
       endif
    end do

    ! Part 1a: BC anthropogenic from GAINS
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(BC_GAINS_filename, BC_GAINS_file, 0)
    call cam_pio_openfile( fid, BC_GAINS_file, PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'time', fid_time)
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_time,f_ntime)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(BC_lat(f_nlat))
    allocate(BC_lon(f_nlon))
    allocate(BC_f3d(f_nlon, f_nlat, f_ntime))
    allocate(BC_f2d(f_nlon, f_nlat))
    allocate(BC_dom_f2d(f_nlon, f_nlat))
    allocate(BC_dom_f3d(f_nlon, f_nlat, f_ntime))
    allocate(BC_awb_f3d(f_nlon, f_nlat, f_ntime))
    allocate(BC2d (plon, plat))
    allocate(BC2d_dom (plon, plat))
    allocate(BC_anthro_GAINS(nmonth, plat, plon))

    ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emis_all', idvar)
    ierr = pio_get_var(fid, idvar, BC_f3d )
    ierr = pio_inq_varid(fid, 'emis_dom', idvar_dom)
    ierr = pio_get_var(fid, idvar, BC_dom_f3d )
    ierr = pio_inq_varid(fid, 'emis_awb', idvar_awb)
    ierr = pio_get_var(fid, idvar, BC_awb_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, BC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, BC_lon )

    ! Close the file.
    call pio_closefile(fid)
    ! get emission excluding domestic and agriculture waste buring
    BC_f2d = BC_f3d(:,:,1) - BC_dom_f3d(:,:,1) - BC_awb_f3d(:,:,1)
    BC_dom_f2d = BC_dom_f3d(:,:,1)

    ! make sure file longitude range from 0-360
    if (BC_lon(1) < -160._r8) then
       allocate(tempor2d(f_nlon, f_nlat))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       ! emission excluding dom
       tempor2d(1:mid_idx,:f_nlat) = BC_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor1d(1:mid_idx) = BC_lon(mid_idx+1:f_nlon)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = BC_f2d(1:mid_idx,:f_nlat)
       tempor1d(mid_idx+1:f_nlon) = BC_lon(1:mid_idx)+360._r8
       BC_f2d = tempor2d
       ! dom emission
       tempor2d(1:mid_idx,:f_nlat) = BC_dom_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = BC_dom_f2d(1:mid_idx,:f_nlat)
       BC_dom_f2d = tempor2d
       !
       BC_lon = tempor1d
       deallocate(tempor2d)
       deallocate(tempor1d)
    else
       BC_lon = BC_lon
    endif

    ! Convert kt/year ----> #/cm2/s
    rearth = 6.371e6_r8 ! m
    do i = 1, f_nlat
       gridarea = 2.0_r8*3.14159_r8*rearth/f_nlat * &
                          2.0_r8*3.14159_r8*rearth/f_nlon*cos(BC_lat(i)/180._r8*3.14159_r8)
       !
       BC_f2d(:f_nlon,i) = BC_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                                        12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
       !
       BC_dom_f2d(:f_nlon,i) = BC_dom_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                                        12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
    end do

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(BC_f2d, f_nlon, f_nlat, BC2d, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(BC_dom_f2d, f_nlon, f_nlat, BC2d_dom, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    ! To implement Monthly data for dom emssion
    ! methods from Stohl et al., 2013
    ! facH works for high latitudes: 45-90N
    ! facL works for low latitudes: 15-45N
    ! below 15N, no seasonal variation
    !
    do itime = 1, nmonth
       ! 45N-90N
       BC2d(:plon, ind_45N:plat) = BC2d(:plon, ind_45N:plat) + &
                                   BC2d_dom(:plon, ind_45N:plat)*facH(itime)*12._r8
       ! 15N-45N
       BC2d(:plon, ind_15N:ind_45N-1) = BC2d(:plon, ind_15N:ind_45N-1) + &
                                        BC2d_dom(:plon, ind_15N:ind_45N-1)*facL(itime)*12._r8
       ! 90S-15N
       BC2d(:plon, 1:ind_15N-1) = BC2d(:plon, 1:ind_15N-1) + &
                                  BC2d_dom(:plon, 1:ind_15N-1)

       BC_anthro_GAINS(itime, :plat, :plon) = transpose(BC2d(:plon, :plat))
    end do

    deallocate(BC_lat)
    deallocate(BC_lon)
    deallocate(BC_f2d)
    deallocate(BC_f3d)
    deallocate(BC_dom_f2d)
    deallocate(BC_dom_f3d)
    deallocate(BC_awb_f3d)
    deallocate(BC2d)
    deallocate(BC2d_dom)

    ! Part 1b: OC anthropogenic from GAINS
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(OC_GAINS_filename, OC_GAINS_file, 0)
    call cam_pio_openfile(fid, trim(OC_GAINS_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'time', fid_time)
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_time,f_ntime)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(OC_lat(f_nlat))
    allocate(OC_lon(f_nlon))
    allocate(OC_f2d(f_nlon, f_nlat))
    allocate(OC_f3d(f_nlon, f_nlat, f_ntime))
    allocate(OC_dom_f2d(f_nlon, f_nlat))
    allocate(OC_dom_f3d(f_nlon, f_nlat, f_ntime))
    allocate(OC_awb_f3d(f_nlon, f_nlat, f_ntime))
    allocate(OC2d (plon, plat))
    allocate(OC2d_dom (plon, plat))
    allocate(OC_anthro_GAINS(nmonth, plat, plon))

    ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emis_all', idvar)
    ierr = pio_get_var(fid, idvar, OC_f3d )
    ierr = pio_inq_varid(fid, 'emis_dom', idvar_dom)
    ierr = pio_get_var(fid, idvar, OC_dom_f3d )
    ierr = pio_inq_varid(fid, 'emis_awb', idvar_awb)
    ierr = pio_get_var(fid, idvar, OC_awb_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, OC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, OC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! get emission excluding domestic and agriculture waste burning
    OC_f2d(:,:) = OC_f3d(:,:,1) - OC_dom_f3d(:,:,1) - OC_awb_f3d(:,:,1)
    OC_dom_f2d = OC_dom_f3d(:,:,1)

    ! make sure file longitude range from -180-180 to 0-360
    if (OC_lon(1) < -160._r8) then
       allocate(tempor2d(f_nlon, f_nlat))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       ! emission excluding dom
       tempor2d(1:mid_idx,:f_nlat) = OC_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor1d(1:mid_idx) = OC_lon(mid_idx+1:f_nlon)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = OC_f2d(1:mid_idx,:f_nlat)
       tempor1d(mid_idx+1:f_nlon) = OC_lon(1:mid_idx)+360._r8
       OC_f2d = tempor2d
       ! dom emission
       tempor2d(1:mid_idx,:f_nlat) = OC_dom_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = OC_dom_f2d(1:mid_idx,:f_nlat)
       OC_dom_f2d = tempor2d
       !
       OC_lon = tempor1d
       deallocate(tempor2d)
       deallocate(tempor1d)
    else
       OC_lon = OC_lon
    endif

    ! Convert kt/year ----> #/cm2/s
    rearth = 6.371e6_r8 ! m
    do i = 1, f_nlat
       gridarea = 2.0_r8*3.14159_r8*rearth/f_nlat * &
                  2.0_r8*3.14159_r8*rearth/f_nlon*cos(OC_lat(i)/180._r8*3.14159_r8)
       !
       OC_f2d(:f_nlon,i) = OC_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                           12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
       !
       OC_dom_f2d(:f_nlon,i) = OC_dom_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                               12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
    end do

    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(OC_f2d, f_nlon, f_nlat, OC2d, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(OC_dom_f2d, f_nlon, f_nlat, OC2d_dom, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    ! To implement Monthly data for dom emssion
    ! methods from Stohl et al., 2013
    ! facH works for high latitudes: 45-90N
    ! facL works for low latitudes: 15-45N
    ! below 15N, no seasonal variation
    !
    do itime = 1, nmonth
       ! 45N-90N
       OC2d(:plon, ind_45N:plat) = OC2d(:plon, ind_45N:plat) + &
                                   OC2d_dom(:plon, ind_45N:plat)*facH(itime)*12._r8
       ! 15N-45N
       OC2d(:plon, ind_15N:ind_45N-1) = OC2d(:plon, ind_15N:ind_45N-1) + &
                                        OC2d_dom(:plon, ind_15N:ind_45N-1)*facL(itime)*12._r8
       ! 90S-15N
       OC2d(:plon, 1:ind_15N-1) = OC2d(:plon, 1:ind_15N-1) + &
                                  OC2d_dom(:plon, 1:ind_15N-1)

       OC_anthro_GAINS(itime, :plat, :plon) = transpose(OC2d(:plon, :plat))
    end do

    deallocate(OC_lat)
    deallocate(OC_lon)
    deallocate(OC_f2d)
    deallocate(OC_f3d)
    deallocate(OC_dom_f2d)
    deallocate(OC_dom_f3d)
    deallocate(OC_awb_f3d)
    deallocate(OC2d)
    deallocate(OC2d_dom)

    ! Part 2a: BC ship
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(BC_ship_filename, BC_ship_file, 0)
    call cam_pio_openfile(fid, trim(BC_ship_file), PIO_NOWRITE)
    !call wrap_open(BC_ship_file, 0, fid)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(BC_lat(f_nlat))
    allocate(BC_lon(f_nlon))
    allocate(BC_f3d(f_nlon, f_nlat, nmonth))
    allocate(BC3d (plon, plat, nmonth))
    allocate(BC_ship_GAINS(nmonth, plat, plon))

   ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emiss_shp', idvar)
    ierr = pio_get_var(fid, idvar, BC_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, BC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, BC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (BC_lon(1) < -160._r8) then
       allocate(tempor3d(f_nlon, f_nlat, nmonth))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = BC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = BC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = BC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = BC_lon(1:mid_idx)+360._r8
       BC_f3d = tempor3d
       BC_lon = tempor1d
       deallocate(tempor3d)
       deallocate(tempor1d)
    else
       BC_lon = BC_lon
    endif

    ! convert unit from kg/m2/s to #/cm2/s
    BC_f3d = BC_f3d*1.e3_r8/1.e4_r8/12._r8*6.02e23_r8

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(BC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       BC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       BC_ship_GAINS(itime, :plat, :plon) = transpose(BC3d(:plon, :plat, itime))
    end do

    deallocate(BC_lat)
    deallocate(BC_lon)
    deallocate(BC_f3d)
    deallocate(BC3d)

    ! Part 2b: OC Ship
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(OC_ship_filename, OC_ship_file, 0)
    call cam_pio_openfile(fid, trim(OC_ship_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(OC_lat(f_nlat))
    allocate(OC_lon(f_nlon))
    allocate(OC_f3d(f_nlon, f_nlat, nmonth))
    allocate(OC3d (plon, plat, nmonth))
    allocate(OC_ship_GAINS(nmonth, plat, plon))

    ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emiss_shp', idvar)
    ierr = pio_get_var(fid, idvar, OC_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, OC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, OC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (OC_lon(1) < -160._r8) then
       allocate(tempor3d(f_nlon, f_nlat, nmonth))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = OC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = OC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = OC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = OC_lon(1:mid_idx)+360._r8
       OC_f3d = tempor3d
       OC_lon = tempor1d
       deallocate(tempor3d)
       deallocate(tempor1d)
    else
       OC_lon = OC_lon
    endif

    ! convert unit from kg/m2/s to #/cm2/s
    OC_f3d = OC_f3d*1.e3_r8/1.e4_r8/12._r8*6.02e23_r8

    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(OC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       OC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       OC_ship_GAINS(itime, :plat, :plon) = transpose(OC3d(:plon, :plat, itime))
    end do

    deallocate(OC_lat)
    deallocate(OC_lon)
    deallocate(OC_f3d)
    deallocate(OC3d)

    ! Part 3a: BC GFEDv3
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(BC_GFEDv3_filename, BC_GFEDv3_file, 0)
    call cam_pio_openfile(fid, trim(BC_GFEDv3_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(BC_lat(f_nlat))
    allocate(BC_lon(f_nlon))
    allocate(BC_f3d(f_nlon, f_nlat, nmonth))
    allocate(tempor3d(f_nlon, f_nlat, nmonth))
    allocate(BC3d (plon, plat, nmonth))
    allocate(BC_GFEDv3(nmonth, plat, plon))

    ! Read in the tables.
    BC_f3d = 0._r8
    ierr = pio_inq_varid(fid, 'emis', idvar)
    ierr = pio_get_var(fid, idvar, tempor3d )
    !call wrap_inq_varid(fid, 'emis', idvar)
    !call wrap_get_var_realx(fid, idvar,  tempor3d)
    BC_f3d = BC_f3d + tempor3d
    ! excluding non-real values
    where (BC_f3d(:,:,:) .ge. 1.e10_r8)
        BC_f3d(:,:,:) = 1.e-30_r8
    end where

    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, BC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, BC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (BC_lon(1) < -160._r8) then
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = BC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = BC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = BC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = BC_lon(1:mid_idx)+360._r8
       BC_f3d = tempor3d
       BC_lon = tempor1d
       deallocate(tempor1d)
    else
       BC_lon = BC_lon
    endif

    ! convert unit from g/m2/month to #/cm2/s
    BC_f3d = BC_f3d/1.e4_r8/30._r8/86400._r8/12._r8*6.02e23_r8

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(BC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       BC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       BC_GFEDv3(itime, :plat, :plon) = transpose(BC3d(:plon, :plat, itime))
    end do

    deallocate(BC_lat)
    deallocate(BC_lon)
    deallocate(BC_f3d)
    deallocate(BC3d)
    deallocate(tempor3d)

    ! Part 3b: OC GFEDv3
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(OC_GFEDv3_filename, OC_GFEDv3_file, 0)
    call cam_pio_openfile(fid, trim(OC_GFEDv3_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    ! write(carma%f_LUNOPRT,*) ''
    ! write(carma%f_LUNOPRT,*) 'f_lon = ', f_nlon
    ! write(carma%f_LUNOPRT,*) 'f_lat = ', f_nlat
    ! write(carma%f_LUNOPRT,*) ''

    allocate(OC_lat(f_nlat))
    allocate(OC_lon(f_nlon))
    allocate(OC_f3d(f_nlon, f_nlat, nmonth))
    allocate(tempor3d(f_nlon, f_nlat, nmonth))
    allocate(OC3d (plon, plat, nmonth))
    allocate(OC_GFEDv3(nmonth, plat, plon))

    ! Read in the tables.
     OC_f3d = 0._r8
    ierr = pio_inq_varid(fid, 'emis', idvar)
    ierr = pio_get_var(fid, idvar, tempor3d )
    !call wrap_inq_varid(fid, 'emis', idvar)
    !call wrap_get_var_realx(fid, idvar,  tempor3d)
    OC_f3d = OC_f3d + tempor3d
    ! excluding non-real values
    where (OC_f3d(:,:,:) .ge. 1.e10_r8)
        OC_f3d(:,:,:) = 1.e-30_r8
    end where

    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, OC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, OC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (OC_lon(1) < -160._r8) then
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = OC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = OC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = OC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = OC_lon(1:mid_idx)+360._r8
       OC_f3d = tempor3d
       OC_lon = tempor1d
       deallocate(tempor1d)
    else
       OC_lon = OC_lon
    endif
    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(OC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       OC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       OC_GFEDv3(itime, :plat, :plon) = transpose(OC3d(:plon, :plat, itime))
    end do

    deallocate(OC_lat)
    deallocate(OC_lon)
    deallocate(OC_f3d)
    deallocate(OC3d)
    deallocate(tempor3d)

! Sum
    do itime = 1, nmonth
       BCnew(:plat, :plon, itime) = BC_anthro_GAINS(itime, :plat, :plon) +  &
             BC_ship_GAINS(itime, :plat, :plon) +  BC_GFEDv3(itime, :plat, :plon)
!
       OCnew(:plat, :plon, itime) = OC_anthro_GAINS(itime, :plat, :plon) +  &
             OC_ship_GAINS(itime, :plat, :plon) +  OC_GFEDv3(itime, :plat, :plon)
    end do
!
    deallocate(BC_anthro_GAINS)
    deallocate(OC_anthro_GAINS)
    deallocate(BC_ship_GAINS)
    deallocate(OC_ship_GAINS)
    deallocate(BC_GFEDv3)
    deallocate(OC_GFEDv3)
    deallocate(facH)
    deallocate(facL)
!
    return
  end subroutine CARMAMODEL_BCOCRead

end module carma_model_mod
