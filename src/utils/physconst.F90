module physconst

! Physical constants.  Use csm_share values whenever available.

   use shr_kind_mod,   only: r8 => shr_kind_r8
   use shr_const_mod,  only: shr_const_g,      shr_const_stebol, shr_const_tkfrz,  &
                             shr_const_mwdair, shr_const_rdair,  shr_const_mwwv,   &
                             shr_const_latice, shr_const_latvap, shr_const_cpdair, &
                             shr_const_rhofw,  shr_const_cpwv,   shr_const_rgas,   &
                             shr_const_karman, shr_const_pstd,   shr_const_rhodair,&
                             shr_const_avogad, shr_const_boltz,  shr_const_cpfw,   &
                             shr_const_rwv,    shr_const_zvir,   shr_const_pi,     &
                             shr_const_rearth, shr_const_sday,   shr_const_cday,   &
                             shr_const_spval,  shr_const_omega,  shr_const_cpvir,  &
                             shr_const_tktrip, shr_const_cpice
   use shr_flux_mod,   only: shr_flux_adjust_constants
   use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
   use cam_abortutils, only: endrun
use constituents,   only: pcnst

implicit none
private
save

public  :: physconst_readnl
public  :: physconst_init
public  :: physconst_update
public  :: physconst_calc_kappav
public  :: composition_init
!
! subroutines to compute thermodynamic quantities
!
! See Lauritzen et al. (2018) for formulaes
!
public  :: get_dp                             ! pressure level thickness from dry dp and dry mixing ratios
public  :: get_pmid_from_dp                   ! full level pressure from dp (approximation depends on dycore)
public  :: get_ps                             ! surface pressure
public  :: get_thermal_energy                 ! thermal energy quantity = dp*cp*T
public  :: get_virtual_temp                   ! virtual temperature
public  :: get_cp                             ! (generalized) heat capacity
public  :: get_cp_dry                         ! (generalized) heat capacity for dry air
public  :: get_sum_species                    ! sum of thermodynamically active species: dp_dry*sum_species=dp
public  :: get_virtual_theta                  ! virtual potential temperature
public  :: get_gz                             ! geopotential
public  :: get_gz_given_dp_Tv_Rdry            ! geopotential (with dp,dry R and Tv as input)
public  :: get_Richardson_number              ! Richardson number at layer interfaces
public  :: get_hydrostatic_static_energy      ! geopotential, dry static energy, and kinetic energy
public  :: get_R_dry                          ! (generalized) dry air gas constant
public  :: get_kappa_dry                      ! (generalized) dry kappa = R_dry/cp_dry
public  :: get_dp_ref                         ! reference pressure layer thickness (include topography)
public  :: get_molecular_diff_coef            ! molecular diffusion and thermal conductivity
public  :: get_molecular_diff_coef_reference  ! reference vertical profile of density, molecular diffusion &
                                              ! and thermal conductivity
public  :: get_rho_dry                        ! dry densisty from temperature (temp) and pressure (dp_dry and tracer)
public  :: get_exner                          ! Exner pressure

! Constants based off share code or defined in physconst

real(r8), public, parameter :: avogad      = shr_const_avogad     ! Avogadro's number (molecules/kmole)
real(r8), public, parameter :: boltz       = shr_const_boltz      ! Boltzman's constant (J/K/molecule)
real(r8), public, parameter :: cday        = shr_const_cday       ! sec in calendar day ~ sec
real(r8), public, parameter :: cpliq       = shr_const_cpfw       ! specific heat of fresh h2o (J/K/kg)
real(r8), public, parameter :: cpice       = shr_const_cpice      ! specific heat of ice (J/K/kg)
real(r8), public, parameter :: karman      = shr_const_karman     ! Von Karman constant
real(r8), public, parameter :: latice      = shr_const_latice     ! Latent heat of fusion (J/kg)
real(r8), public, parameter :: latvap      = shr_const_latvap     ! Latent heat of vaporization (J/kg)
real(r8), public, parameter :: pi          = shr_const_pi         ! 3.14...
#ifdef planet_mars
real(r8), public, parameter :: pstd        = 6.0E1_r8             ! Standard pressure (Pascals)
#else
real(r8), public, parameter :: pstd        = shr_const_pstd       ! Standard pressure (Pascals)
real(r8), public, parameter :: tref        = 288._r8              ! Reference temperature
real(r8), public, parameter :: lapse_rate  = 0.0065_r8            ! reference lapse rate [K/m]
#endif
real(r8), public, parameter :: r_universal = shr_const_rgas       ! Universal gas constant (J/K/kmol)
real(r8), public, parameter :: rhoh2o      = shr_const_rhofw      ! Density of liquid water (STP)
real(r8), public, parameter :: spval       = shr_const_spval      !special value
real(r8), public, parameter :: stebol      = shr_const_stebol     ! Stefan-Boltzmann's constant (W/m^2/K^4)
real(r8), public, parameter :: h2otrip     = shr_const_tktrip     ! Triple point temperature of water (K)

real(r8), public, parameter :: c0          = 2.99792458e8_r8      ! Speed of light in a vacuum (m/s)
real(r8), public, parameter :: planck      = 6.6260755e-34_r8     ! Planck's constant (J.s)

! Molecular weights
real(r8), public, parameter :: mwco2       =  44._r8             ! molecular weight co2
real(r8), public, parameter :: mwn2o       =  44._r8             ! molecular weight n2o
real(r8), public, parameter :: mwch4       =  16._r8             ! molecular weight ch4
real(r8), public, parameter :: mwf11       = 136._r8             ! molecular weight cfc11
real(r8), public, parameter :: mwf12       = 120._r8             ! molecular weight cfc12
real(r8), public, parameter :: mwo3        =  48._r8             ! molecular weight O3
real(r8), public, parameter :: mwso2       =  64._r8
real(r8), public, parameter :: mwso4       =  96._r8
real(r8), public, parameter :: mwh2o2      =  34._r8
real(r8), public, parameter :: mwdms       =  62._r8
real(r8), public, parameter :: mwnh4       =  18._r8


! modifiable physical constants for aquaplanet

real(r8), public, protected :: gravit  = shr_const_g      ! gravitational acceleration (m/s**2)
real(r8), public, protected :: sday    = shr_const_sday   ! sec in siderial day ~ sec
real(r8), public, protected :: mwh2o   = shr_const_mwwv   ! molecular weight h2o
real(r8), public, protected :: cpwv    = shr_const_cpwv   ! specific heat of water vapor (J/K/kg)
real(r8), public, protected :: mwdry   = shr_const_mwdair ! molecular weight dry air
real(r8), public, protected :: cpair   = shr_const_cpdair ! specific heat of dry air (J/K/kg)
real(r8), public, protected :: rearth  = shr_const_rearth ! radius of earth (m)
real(r8), public, protected :: tmelt   = shr_const_tkfrz  ! Freezing point of water (K)

!---------------  Variables below here are derived from those above -----------------------

real(r8), public, protected :: rga        = 1._r8/shr_const_g      ! reciprocal of gravit
real(r8), public, protected :: ra         = 1._r8/shr_const_rearth ! reciprocal of earth radius
real(r8), public, protected :: omega      = shr_const_omega        ! earth rot ~ rad/sec
real(r8), public, protected :: rh2o       = shr_const_rwv          ! Water vapor gas constant ~ J/K/kg
real(r8), public, protected :: rair       = shr_const_rdair        ! Dry air gas constant     ~ J/K/kg
real(r8), public, protected :: epsilo     = shr_const_mwwv/shr_const_mwdair   ! ratio of h2o to dry air molecular weights
real(r8), public, protected :: zvir       = shr_const_zvir         ! (rh2o/rair) - 1
real(r8), public, protected :: cpvir      = shr_const_cpvir        ! CPWV/CPDAIR - 1.0
real(r8), public, protected :: rhodair    = shr_const_rhodair      ! density of dry air at STP  ~ kg/m^3
real(r8), public, protected :: cappa      = (shr_const_rgas/shr_const_mwdair)/shr_const_cpdair  ! R/Cp
real(r8), public, protected :: ez         ! Coriolis expansion coeff -> omega/sqrt(0.375)
real(r8), public, protected :: Cpd_on_Cpv = shr_const_cpdair/shr_const_cpwv

!---------------  Variables below here are for WACCM-X -----------------------
real(r8), public, dimension(:,:,:), pointer :: cpairv ! composition dependent specific heat at constant pressure
real(r8), public, dimension(:,:,:), pointer :: rairv  ! composition dependent gas "constant"
real(r8), public, dimension(:,:,:), pointer :: cappav ! rairv/cpairv
real(r8), public, dimension(:,:,:), pointer :: mbarv  ! composition dependent atmosphere mean mass
real(r8), public, dimension(:,:,:), pointer :: kmvis  ! molecular viscosity      kg/m/s
real(r8), public, dimension(:,:,:), pointer :: kmcnd  ! molecular conductivity   J/m/s/K

!--------------- Variables for consistent themodynamics --------------------
!
! composition of air
!
integer, parameter :: num_names_max = 30
character(len=6 )  :: dry_air_species(num_names_max)
character(len=6 )  :: water_species_in_air(num_names_max)

integer, protected, public      :: dry_air_species_num
integer, protected, public      :: water_species_in_air_num

integer,               protected, public :: thermodynamic_active_species_num
integer,  allocatable, protected, public :: thermodynamic_active_species_idx(:)
integer,  allocatable,            public :: thermodynamic_active_species_idx_dycore(:)    
real(r8), allocatable, protected, public :: thermodynamic_active_species_cp(:)
real(r8), allocatable, protected, public :: thermodynamic_active_species_cv(:)
real(r8), allocatable, protected, public :: thermodynamic_active_species_R(:)
real(r8), allocatable, protected, public :: thermodynamic_active_species_mwi(:)!inverse molecular weights dry air
real(r8), allocatable, protected, public :: thermodynamic_active_species_kv(:) !molecular diffusion
real(r8), allocatable, protected, public :: thermodynamic_active_species_kc(:) !thermal conductivity

! standard dry air (constant composition)
real(r8) :: mmro2, mmrn2           ! Mass mixing ratios of O2 and N2
real(r8) :: o2_mwi, n2_mwi         ! Inverse molecular weights
real(r8) :: mbar                   ! Mean mass at mid level

! coefficients in expressions for molecular diffusion coefficients
! kv1,..,kv4 are coefficients for kmvis calculation
! kc1,..,kc4 are coefficients for kmcnd calculation
real(r8), parameter :: &
   kv1 = 4.03_r8, &
   kv2 = 3.42_r8, &
   kv3 = 3.9_r8,  &
   kv4 = 0.69_r8, &
   kc1 = 56._r8,  &
   kc2 = 56._r8,  &
   kc3 = 75.9_r8, &
   kc4 = 0.69_r8

!================================================================================================
contains
!================================================================================================

! Read namelist variables.
subroutine physconst_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use spmd_utils,      only: masterproc, mpicom, masterprocid, mpi_real8, mpi_character
   use cam_logfile,     only: iulog

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=*), parameter :: subname = 'physconst_readnl'
   logical :: newg, newsday, newmwh2o, newcpwv, newmwdry, newcpair, newrearth, newtmelt, newomega


   ! Physical constants needing to be reset (e.g., for aqua planet experiments)
   namelist /physconst_nl/  gravit, sday, mwh2o, cpwv, mwdry, cpair, rearth, tmelt, omega

   ! Variable components of dry air and water species in air
   namelist /air_composition_nl/ dry_air_species, water_species_in_air
   !-----------------------------------------------------------------------------

   if (masterproc) then
      open( newunit=unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'physconst_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, physconst_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
   end if

   ! Broadcast namelist variables
   call MPI_bcast(gravit, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(sday,   1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(mwh2o,  1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(cpwv,   1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(mwdry,  1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(cpair,  1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(rearth, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(tmelt,  1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(omega,  1, mpi_real8, masterprocid, mpicom, ierr)

   newg     =  gravit .ne. shr_const_g
   newsday  =  sday   .ne. shr_const_sday
   newmwh2o =  mwh2o  .ne. shr_const_mwwv
   newcpwv  =  cpwv   .ne. shr_const_cpwv
   newmwdry =  mwdry  .ne. shr_const_mwdair
   newcpair =  cpair  .ne. shr_const_cpdair
   newrearth=  rearth .ne. shr_const_rearth
   newtmelt =  tmelt  .ne. shr_const_tkfrz
   newomega =  omega  .ne. shr_const_omega

   if (newg .or. newsday .or. newmwh2o .or. newcpwv .or. newmwdry .or. newrearth .or. newtmelt .or. newomega) then
      if (masterproc) then
         write(iulog,*)'****************************************************************************'
         write(iulog,*)'***    New Physical Constant Values set via namelist                     ***'
         write(iulog,*)'***                                                                      ***'
         write(iulog,*)'***    Physical Constant    Old Value                  New Value         ***'
         if (newg)       write(iulog,*)'***       GRAVIT    ',shr_const_g,gravit,'***'
         if (newsday)    write(iulog,*)'***       SDAY      ',shr_const_sday,sday,'***'
         if (newmwh2o)   write(iulog,*)'***       MWH20     ',shr_const_mwwv,mwh2o,'***'
         if (newcpwv)    write(iulog,*)'***       CPWV      ',shr_const_cpwv,cpwv,'***'
         if (newmwdry)   write(iulog,*)'***       MWDRY     ',shr_const_mwdair,mwdry,'***'
         if (newcpair)   write(iulog,*)'***       CPAIR     ',shr_const_cpdair,cpair,'***'
         if (newrearth)  write(iulog,*)'***       REARTH    ',shr_const_rearth,rearth,'***'
         if (newtmelt)   write(iulog,*)'***       TMELT     ',shr_const_tkfrz,tmelt,'***'
         if (newomega)   write(iulog,*)'***       OMEGA     ',shr_const_omega,omega,'***'
         write(iulog,*)'****************************************************************************'
      end if
      rga         = 1._r8/gravit
      ra          = 1._r8/rearth
      if (.not. newomega) then
         omega       = 2.0_r8*pi/sday
      end if
      cpvir       = cpwv/cpair - 1._r8
      epsilo      = mwh2o/mwdry

      !  rair and rh2o have to be defined before any of the variables that use them

      rair        = r_universal/mwdry
      rh2o        = r_universal/mwh2o

      cappa       = rair/cpair
      rhodair     = pstd/(rair*tmelt)
      zvir        =  (rh2o/rair)-1.0_R8
      ez          = omega / sqrt(0.375_r8)
      Cpd_on_Cpv  = cpair/cpwv

      ! Adjust constants in shr_flux_mod.
      call shr_flux_adjust_constants(zvir=zvir, cpvir=cpvir, gravit=gravit)

   else
      ez          = omega / sqrt(0.375_r8)
   end if

   ! Read variable components of dry air and water species in air

   dry_air_species   = (/ (' ', i=1,num_names_max) /)
   water_species_in_air = (/ (' ', i=1,num_names_max) /)    
    
   if (masterproc) then
      open( newunit=unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'air_composition_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, air_composition_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
   end if
    
   call mpi_bcast(dry_air_species, len(dry_air_species)*num_names_max, mpi_character, &
                  masterprocid, mpicom, ierr)
   call mpi_bcast(water_species_in_air, len(water_species_in_air)*num_names_max, mpi_character, &
                  masterprocid, mpicom, ierr)
    
   dry_air_species_num   = 0
   water_species_in_air_num = 0
   do i = 1, num_names_max
      if (.not. LEN(TRIM(dry_air_species(i)))==0) then
         dry_air_species_num = dry_air_species_num + 1
      end if
      if (.not. LEN(TRIM(water_species_in_air(i)))==0) then
         water_species_in_air_num = water_species_in_air_num + 1
      endif
   end do
   thermodynamic_active_species_num = dry_air_species_num+water_species_in_air_num

   if (masterproc) then

      write(iulog,*)'****************************************************************************'
      write(iulog,*)' '

      if (dry_air_species_num == 0) then
         write(iulog,*)' Thermodynamic properties of dry air are fixed at troposphere values'
      else
         write(iulog,*)' Thermodynamic properties of dry air are based on variable'
         write(iulog,*)' composition of the following species:'
         do i = 1, dry_air_species_num
            write(iulog,*)'   ', trim(dry_air_species(i))
         end do
         write(iulog,*) ' '
      end if

      write(iulog,*)' Thermodynamic properties of moist air are based on variable'
      write(iulog,*)' composition of the following water species:'
      do i = 1, water_species_in_air_num
         write(iulog,*)'   ', trim(water_species_in_air(i))
      end do

      write(iulog,*)' '
      write(iulog,*)'****************************************************************************'

   end if

end subroutine physconst_readnl

!===============================================================================

subroutine physconst_init()

   integer :: ierr

   !-------------------------------------------------------------------------------
   !  Allocate constituent dependent properties
   !-------------------------------------------------------------------------------
   allocate( cpairv(pcols,pver,begchunk:endchunk), &
             rairv(pcols,pver,begchunk:endchunk),  &
             cappav(pcols,pver,begchunk:endchunk), &
             mbarv(pcols,pver,begchunk:endchunk),  &
             kmvis(pcols,pverp,begchunk:endchunk), &
             kmcnd(pcols,pverp,begchunk:endchunk), stat=ierr )
   if ( ierr /= 0 ) call endrun('physconst: allocate failed in physconst_init')

   !-------------------------------------------------------------------------------
   !  Initialize constituent dependent properties
   !-------------------------------------------------------------------------------
   cpairv(:pcols,:pver,begchunk:endchunk) = cpair
   rairv(:pcols,:pver,begchunk:endchunk) = rair
   cappav(:pcols,:pver,begchunk:endchunk) = rair/cpair
   mbarv(:pcols,:pver,begchunk:endchunk) = mwdry

end subroutine physconst_init

!===============================================================================

  subroutine composition_init()
    use constituents, only: cnst_get_ind, cnst_mw
    use spmd_utils,      only: masterproc
    use cam_logfile,     only: iulog
    character(len=*), parameter :: subname = 'composition_init'
    real(r8) :: mw, dof1, dof2, dof3

    integer  :: icnst,ix,i

    ! standard dry air (constant composition)
    o2_mwi = 1._r8/32._r8
    n2_mwi = 1._r8/28._r8
    mmro2  = 0.235_r8
    mmrn2  = 0.765_r8
    mbar = 1._r8/(mmro2*o2_mwi + mmrn2*n2_mwi)
    
    ! init for variable composition dry air

    i = dry_air_species_num+water_species_in_air_num  
    allocate(thermodynamic_active_species_idx(i))
    allocate(thermodynamic_active_species_idx_dycore(i))
    allocate(thermodynamic_active_species_cp(0:i))
    allocate(thermodynamic_active_species_cv(0:i))
    allocate(thermodynamic_active_species_R(0:i))
    
    i = dry_air_species_num
    allocate(thermodynamic_active_species_mwi(i))
    allocate(thermodynamic_active_species_kv(i))
    allocate(thermodynamic_active_species_kc(i))
    thermodynamic_active_species_idx        = -999
    thermodynamic_active_species_idx_dycore = -999
    thermodynamic_active_species_cp         = 0.0_r8
    thermodynamic_active_species_cv         = 0.0_r8
    thermodynamic_active_species_R          = 0.0_r8
    thermodynamic_active_species_mwi        = 0.0_r8
    thermodynamic_active_species_kv         = 0.0_r8
    thermodynamic_active_species_kc         = 0.0_r8
    !
    ! define cp and R for species in species_name
    !
    ! Last major species in namelist dry_air_species is derived from the other major species
    ! (since sum of dry mixing ratios for major species of dry air add must add to one)
    !
    dof1 = 3._r8   ! monatomic ideal gas  cv=dof1/2 * R;   cp=(1+dof1/2) * R; dof=3 translational
    dof2 = 5._r8   ! diatomic ideal gas   cv=dof2/2 * R;   cp=(1+dof2/2) * R; dof=3 tranlational + 2 rotational
    dof3 = 6._r8   ! polyatomic ideal gas cv=dof3/2 * R;   cp=(1+dof3/2) * R; dof=3 tranlational + 3 rotational
    !
    if (dry_air_species_num>0) then
      !
      ! last major species in dry_air_species is derived from the others and constants associated with it
      ! are initialized here
      !
      if (TRIM(dry_air_species(dry_air_species_num))=='N2') then
        call cnst_get_ind('N' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' dry air component not found: ', dry_air_species(dry_air_species_num)
          call endrun(subname // ':: dry air component not found')         
        else         
          mw = 2.0_r8*cnst_mw(ix)
          icnst = dry_air_species_num
          thermodynamic_active_species_idx(icnst) = 1!note - this is not used since this tracer value is derived
          thermodynamic_active_species_cp (icnst) = 0.5_r8*shr_const_rgas*(2._r8+dof2)/mw !N2
          thermodynamic_active_species_cv (icnst) = 0.5_r8*shr_const_rgas*dof2/mw !N2
          thermodynamic_active_species_R  (icnst) = shr_const_rgas/mw
          thermodynamic_active_species_mwi(icnst) = 1.0_r8/mw
          thermodynamic_active_species_kv(icnst)  = 3.42_r8
          thermodynamic_active_species_kc(icnst)  = 56._r8
        end if
        !
        ! if last major species is not N2 then add code here
        !
      else
        write(iulog, *) subname//' derived major species not found: ', dry_air_species(dry_air_species_num)
        call endrun(subname // ':: derived major species not found')
      end if
    else
      !
      ! dry air is not species dependent
      !
      icnst = 0
      thermodynamic_active_species_cp (icnst) = cpair
      thermodynamic_active_species_cv (icnst) = cpair - rair
      thermodynamic_active_species_R  (icnst) = rair
    end if
    !
    !******************************************************************************
    !
    ! add prognostic components of dry air
    !
    !******************************************************************************
    !
    icnst = 1
    do i=1,dry_air_species_num-1
      select case (TRIM(dry_air_species(i)))
        !
        ! O
        !
      case('O')
        call cnst_get_ind('O' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' dry air component not found: ', dry_air_species(i)
          call endrun(subname // ':: dry air component not found')         
        else
          mw = cnst_mw(ix)
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = 0.5_r8*shr_const_rgas*(2._r8+dof1)/mw
          thermodynamic_active_species_cv (icnst) = 0.5_r8*shr_const_rgas*dof1/mw
          thermodynamic_active_species_R  (icnst) = shr_const_rgas/mw
          thermodynamic_active_species_mwi(icnst) = 1.0_r8/mw
          thermodynamic_active_species_kv(icnst)  = 3.9_r8
          thermodynamic_active_species_kc(icnst)  = 75.9_r8
          icnst = icnst+1
        end if
        !
        ! O2
        !
      case('O2')
        call cnst_get_ind('O2' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' dry air component not found: ', dry_air_species(i)
          call endrun(subname // ':: dry air component not found')         
        else
          mw = cnst_mw(ix)
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = 0.5_r8*shr_const_rgas*(2._r8+dof2)/mw
          thermodynamic_active_species_cv (icnst) = 0.5_r8*shr_const_rgas*dof2/mw
          thermodynamic_active_species_R  (icnst) = shr_const_rgas/mw
          thermodynamic_active_species_mwi(icnst) = 1.0_r8/mw
          thermodynamic_active_species_kv(icnst)  = 4.03_r8
          thermodynamic_active_species_kc(icnst)  = 56._r8
          icnst = icnst+1
        end if
        !
        ! H
        !
      case('H')
        call cnst_get_ind('H' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' dry air component not found: ', dry_air_species(i)
          call endrun(subname // ':: dry air component not found')         
        else
          mw = cnst_mw(ix)
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = 0.5_r8*shr_const_rgas*(2._r8+dof1)/mw
          thermodynamic_active_species_cv (icnst) = 0.5_r8*shr_const_rgas*dof1/mw
          thermodynamic_active_species_R  (icnst) = shr_const_rgas/mw
          thermodynamic_active_species_mwi(icnst) = 1.0_r8/mw
          thermodynamic_active_species_kv(icnst)  = 0.0_r8
          thermodynamic_active_species_kc(icnst)  = 0.0_r8
          icnst = icnst+1
        end if
        !
        ! If support for more major species is to be included add code here
        !       
      case default
        write(iulog, *) subname//' dry air component not found: ', dry_air_species(i)
        call endrun(subname // ':: dry air component not found')         
      end select
      
      if (masterproc) then
        write(iulog, *) "Dry air composition ",TRIM(dry_air_species(i)),&
             icnst-1,thermodynamic_active_species_idx(icnst-1),&
             thermodynamic_active_species_mwi(icnst-1),&
             thermodynamic_active_species_cp(icnst-1),&
             thermodynamic_active_species_cv(icnst-1)
      end if
    end do
    i = dry_air_species_num
    if (i>0) then
      if (masterproc) then
        write(iulog, *) "Dry air composition ",TRIM(dry_air_species(i)),&
             icnst,thermodynamic_active_species_idx(icnst),&
             thermodynamic_active_species_mwi(icnst),&
             thermodynamic_active_species_cp(icnst),&
             thermodynamic_active_species_cv(icnst)
      end if
    end if
    !
    !************************************************************************************
    !
    ! Add non-dry components of moist air (water vapor and condensates)
    !
    !************************************************************************************
    !
    icnst = dry_air_species_num+1
    do i=1,water_species_in_air_num
      select case (TRIM(water_species_in_air(i)))
        !
        ! Q
        !
      case('Q')
        call cnst_get_ind('Q' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' moist air component not found: ', water_species_in_air(i)
          call endrun(subname // ':: moist air component not found')         
        else
          mw = cnst_mw(ix)
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = cpwv
          thermodynamic_active_species_cv (icnst) = 0.5_r8*shr_const_rgas*dof3/mw
          thermodynamic_active_species_R  (icnst) = rh2o
          icnst = icnst+1
        end if
        !
        ! CLDLIQ
        !
      case('CLDLIQ')
        call cnst_get_ind('CLDLIQ' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' moist air component not found: ', water_species_in_air(i)
          call endrun(subname // ':: moist air component not found')         
        else
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = cpliq
          thermodynamic_active_species_cv (icnst) = cpliq
          icnst = icnst+1
        end if
        !
        ! CLDICE
        !
      case('CLDICE')
        call cnst_get_ind('CLDICE' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' moist air component not found: ', water_species_in_air(i)
          call endrun(subname // ':: moist air component not found')         
        else
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = cpice
          thermodynamic_active_species_cv (icnst) = cpice
          icnst = icnst+1
        end if
        !
        ! RAINQM
        !
      case('RAINQM')
        call cnst_get_ind('RAINQM' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' moist air component not found: ', water_species_in_air(i)
          call endrun(subname // ':: moist air component not found')         
        else
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = cpliq
          thermodynamic_active_species_cv (icnst) = cpliq
          icnst = icnst+1
        end if
        !
        ! SNOWQM
        !
      case('SNOWQM')
        call cnst_get_ind('SNOWQM' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' moist air component not found: ', water_species_in_air(i)
          call endrun(subname // ':: moist air component not found')         
        else
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = cpice
          thermodynamic_active_species_cv (icnst) = cpice
          icnst = icnst+1
        end if
        !
        ! GRAUQM
        !
      case('GRAUQM')
        call cnst_get_ind('GRAUQM' ,ix, abort=.false.)
        if (ix<1) then
          write(iulog, *) subname//' moist air component not found: ', water_species_in_air(i)
          call endrun(subname // ':: moist air component not found')         
        else
          mw = cnst_mw(ix)
          thermodynamic_active_species_idx(icnst) = ix
          thermodynamic_active_species_cp (icnst) = cpice
          thermodynamic_active_species_cv (icnst) = cpice
          icnst = icnst+1
        end if
        !
        ! If support for more major species is to be included add code here
        !       
      case default
        write(iulog, *) subname//' moist air component not found: ', water_species_in_air(i)
        call endrun(subname // ':: moist air component not found')         
      end select
      !
      !
      !
      if (masterproc) then
        write(iulog, *) "Thermodynamic active species ",TRIM(water_species_in_air(i)),&
             icnst-1,thermodynamic_active_species_idx(icnst-1),&
             thermodynamic_active_species_cp(icnst-1),&
             thermodynamic_active_species_cv(icnst-1)
      end if
    end do
    
  end subroutine composition_init  
  !
  !****************************************************************************************************************
  !
  ! update species dependent constants for physics
  !
  !****************************************************************************************************************
  ! 
  subroutine physconst_update(mmr, t, lchnk, ncol, to_moist_factor)
    
    !-----------------------------------------------------------------------
    ! Update the physics "constants" that vary
    !-----------------------------------------------------------------------
    
    !------------------------------Arguments--------------------------------------------------------------
    
    real(r8), intent(in) :: mmr(pcols,pver,pcnst) ! constituents q array from state structure
    real(r8), intent(in) :: t(pcols,pver)   ! temperature t array from state structure
    integer, intent(in)  :: lchnk           ! Chunk number
    integer, intent(in)  :: ncol            ! number of columns
    real(r8),  optional, intent(in) :: to_moist_factor(:,:)
    !
    !---------------------------Local storage-------------------------------------------------------------
    real(r8):: to_moist_fact(ncol,pver)
    real(r8):: sponge_factor(pver)
    
    to_moist_fact(:,:) = 1._r8
    
    if (present(to_moist_factor)) then
      to_moist_fact(:ncol,:) = to_moist_factor(:ncol,:)
    end if
    
    !--------------------------------------------
    ! update cpairv, rairv, mbarv, and cappav
    !--------------------------------------------
    call get_R_dry(1,ncol,1,1,1,pver,1,pver,pcnst, mmr(:ncol,:,:), thermodynamic_active_species_idx, &
                   rairv(:ncol,:,lchnk), fact=to_moist_fact(:ncol,:))
    call get_cp_dry(1,ncol,1,1,1,pver,1,pver,pcnst, mmr(:ncol,:,:), thermodynamic_active_species_idx, &
                    cpairv(:ncol,:,lchnk), fact=to_moist_fact(:ncol,:))
    call get_mbarv(1,ncol,1,1,1,pver,pver,pcnst, mmr(:ncol,:,:), thermodynamic_active_species_idx, &
                   mbarv(:ncol,:,lchnk), fact=to_moist_fact(:ncol,:))
    sponge_factor = 1.0_r8
    call get_molecular_diff_coef(1,ncol,1,1,pver,pver,t(:ncol,:),1,sponge_factor,kmvis(:ncol,:,lchnk), &
       kmcnd(:ncol,:,lchnk), pcnst, tracer=mmr(:ncol,:,:), fact=to_moist_fact(:ncol,:),                &
       active_species_idx_dycore=thermodynamic_active_species_idx)

  end subroutine physconst_update
  !
  !****************************************************************************************************************
  !
  ! update species dependent kappa for FV dycore
  !
  !****************************************************************************************************************
  ! 
  subroutine physconst_calc_kappav( i0,i1,j0,j1,k0,k1,ntotq, tracer, kappav, cpv )
    ! assumes moist MMRs
    
    ! args
    integer,  intent(in) :: i0,i1,j0,j1,k0,k1, ntotq
    real(r8), intent(in) :: tracer(i0:i1,j0:j1,k0:k1,ntotq) ! Tracer array
     real(r8), intent(out) :: kappav(i0:i1,j0:j1,k0:k1)
     real(r8), optional, intent(out) :: cpv(i0:i1,j0:j1,k0:k1) 

     ! local vars
     real(r8),  dimension(i0:i1,j0:j1,k0:k1) :: rgas_var, cp_var
     integer :: i,j,k

     !-----------------------------------------------------------------------
     !  Calculate constituent dependent specific heat, gas constant and cappa
     !-----------------------------------------------------------------------
     call get_R_dry (i0,i1,j0,j1,k0,k1,k0,k1,ntotq,tracer,thermodynamic_active_species_idx,rgas_var)
     call get_cp_dry(i0,i1,j0,j1,k0,k1,k0,k1,ntotq,tracer,thermodynamic_active_species_idx,cp_var)
!$omp parallel do private(i,j,k)
     do k = k0,k1
        do j = j0,j1
          do i = i0,i1
            kappav(i,j,k)   = rgas_var(i,j,k)/cp_var(i,j,k)
           enddo
        enddo
     enddo

     if (present(cpv)) then
        cpv(:,:,:) = cp_var(:,:,:)
     endif

   end subroutine physconst_calc_kappav
   !
   !****************************************************************************************************************
   !
   ! Compute pressure level thickness from dry pressure and thermodynamic active species mixing ratios
   !
   ! Tracer can either be in units of dry mixing ratio (mixing_ratio=1) or "mass" (=m*dp_dry) (mixing_ratio=2)
   !
   !****************************************************************************************************************
   !
   subroutine get_dp(i0,i1,j0,j1,k0,k1,ntrac,tracer,mixing_ratio,active_species_idx,dp_dry,dp,ps,ptop)
     integer,  intent(in)  :: i0,i1,j0,j1,k0,k1,ntrac            ! array bounds
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,k0:k1,1:ntrac)  ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is dry mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array     
     real(r8), intent(in)  :: dp_dry(i0:i1,j0:j1,k0:k1)          ! dry pressure level thickness
     real(r8), intent(out) :: dp(i0:i1,j0:j1,k0:k1)              ! pressure level thickness
     real(r8), optional,intent(out) :: ps(i0:i1,j0:j1)           ! surface pressure (if ps present then ptop
                                                                 !                   must be present)
     real(r8), optional,intent(in)  :: ptop                      ! pressure at model top

     integer :: i,j,k,m_cnst,nq
    
     dp = dp_dry
     if (mixing_ratio==1) then
       do nq=dry_air_species_num+1,thermodynamic_active_species_num       
         m_cnst = active_species_idx(nq)       
         do k=k0,k1
           do j=j0,j1
             do i = i0,i1
               dp(i,j,k) = dp(i,j,k) + dp_dry(i,j,k)*tracer(i,j,k,m_cnst)
             end do
           end do
         end do
       end do       
     else
       do nq=dry_air_species_num+1,thermodynamic_active_species_num              
         m_cnst = active_species_idx(nq)       
         do k=k0,k1
           do j=j0,j1
             do i = i0,i1
               dp(i,j,k) = dp(i,j,k) + tracer(i,j,k,m_cnst)
             end do
           end do
         end do
       end do
     end if
     if (present(ps)) then
       if (present(ptop)) then
         ps = ptop
         do k=k0,k1
           do j=j0,j1
             do i = i0,i1
               ps(i,j) = ps(i,j)+dp(i,j,k)
             end do
           end do
         end do
       else
         call endrun('get_dp: if ps is present ptop must be present')
       end if
     end if
   end subroutine get_dp
   !
   !*************************************************************************************************************************
   !
   ! compute mid-level (full level) pressure from dry pressure and water tracers
   !
   !*************************************************************************************************************************
   !
   subroutine get_pmid_from_dpdry(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                                  dp_dry, ptop, pmid, pint, dp)

     integer,  intent(in)  :: i0,i1,j0,j1,nlev,ntrac             ! array bounds
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)   ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(i0:i1,j0:j1,nlev)           ! dry pressure level thickness     
     real(r8), intent(in)  :: ptop                               ! model top pressure
     real(r8), intent(out) :: pmid(i0:i1,j0:j1,nlev)             ! mid-level pressure
     real(r8), optional, intent(out) :: pint(i0:i1,j0:j1,nlev+1) ! half-level pressure
     real(r8), optional, intent(out) :: dp(i0:i1,j0:j1,nlev)     ! presure level thickness

     real(r8) :: dp_local(i0:i1,j0:j1,nlev)                      ! local pressure level thickness
     real(r8) :: pint_local(i0:i1,j0:j1,nlev+1)                  ! local interface pressure
     integer  :: k

     call get_dp(i0,i1,j0,j1,1,nlev,ntrac,tracer,mixing_ratio,active_species_idx,dp_dry,dp_local)
     pint_local(:,:,1) = ptop
     do k=2,nlev+1
       pint_local(:,:,k) = dp_local(:,:,k-1)+pint_local(:,:,k-1)
     end do

     call get_pmid_from_dp(i0,i1,j0,j1,1,nlev,dp_local,ptop,pmid,pint_local)
     
     if (present(pint)) pint=pint_local
     if (present(dp)) dp=dp_local
   end subroutine get_pmid_from_dpdry
   !
   !*************************************************************************************************************************
   !
   ! compute mid-level (full level) pressure
   !
   !*************************************************************************************************************************
   !
   subroutine get_pmid_from_dp(i0,i1,j0,j1,k0,k1,dp,ptop,pmid,pint)
     use dycore, only: dycore_is     
     integer,  intent(in)            :: i0,i1,j0,j1,k0,k1         ! array bounds
     real(r8), intent(in)            :: dp(i0:i1,j0:j1,k0:k1)     ! dry pressure level thickness     
     real(r8), intent(in)            :: ptop                      ! pressure at model top
     real(r8), intent(out)           :: pmid(i0:i1,j0:j1,k0:k1)   ! mid (full) level pressure
     real(r8), optional, intent(out) :: pint(i0:i1,j0:j1,k0:k1+1) ! pressure at interfaces (half levels)

     real(r8) :: pint_local(i0:i1,j0:j1,k0:k1+1)        
     integer  :: k

     pint_local(:,:,k0) = ptop
     do k=k0+1,k1+1
       pint_local(:,:,k) = dp(:,:,k-1)+pint_local(:,:,k-1)
     end do

     if (dycore_is ('LR')) then
       do k=k0,k1
         pmid(:,:,k) = dp(:,:,k)/(log(pint_local(:,:,k+1))-log(pint_local(:,:,k)))
       end do
     else
       do k=k0,k1
         pmid(:,:,k) = 0.5_r8*(pint_local(:,:,k)+pint_local(:,:,k+1))
       end do       
     end if
     if (present(pint)) pint=pint_local     
   end subroutine get_pmid_from_dp   
   !
   !****************************************************************************************************************   
   !   
   ! Compute Exner pressure
   !
   !****************************************************************************************************************   
   !
   subroutine get_exner(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx,&
        dp_dry,ptop,p00,inv_exner,exner,poverp0)
     integer,  intent(in)  :: i0,i1,j0,j1,nlev,ntrac             ! index bounds
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)   ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array          
     real(r8), intent(in)  :: dp_dry(i0:i1,j0:j1,nlev)           ! dry pressure level thickness     
     real(r8), intent(in)  :: ptop                               ! pressure at model top
     real(r8), intent(in)  :: p00                                ! reference pressure for Exner pressure (usually 1000hPa)
     logical , intent(in)  :: inv_exner                          ! logical for outputting inverse Exner or Exner pressure
     real(r8), intent(out) :: exner(i0:i1,j0:j1,nlev)        
     real(r8), optional, intent(out) :: poverp0(i0:i1,j0:j1,nlev)! for efficiency when a routine needs this variable

     real(r8) :: pmid(i0:i1,j0:j1,nlev),kappa_dry(i0:i1,j0:j1,nlev)
     !
     ! compute mid level pressure
     !
     call get_pmid_from_dpdry(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx,dp_dry,ptop,pmid)
     !
     ! compute kappa = Rd/cpd
     !
     if (mixing_ratio==1) then
       call get_kappa_dry(i0,i1,j0,j1,1,nlev,nlev,ntrac,tracer,active_species_idx,kappa_dry)
     else
       call get_kappa_dry(i0,i1,j0,j1,1,nlev,nlev,ntrac,tracer,active_species_idx,kappa_dry,1.0_r8/dp_dry)
     end if
     if (inv_exner) then
       exner(:,:,:) = (p00/pmid(:,:,:))**kappa_dry(:,:,:)       
     else
       exner(:,:,:) = (pmid(:,:,:)/p00)**kappa_dry(:,:,:)
     end if
     if (present(poverp0)) poverp0=pmid(:,:,:)/p00
   end subroutine get_exner
   !
   !****************************************************************************************************************   
   !   
   ! Compute virtual potential temperature from dp_dry, m, T and ptop.
   !
   !****************************************************************************************************************   
   !   
   subroutine get_virtual_theta(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                                dp_dry,ptop,p00,temp,theta_v)
     integer,  intent(in)  :: i0,i1,j0,j1,nlev,ntrac
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)   ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is dry mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(i0:i1,j0:j1,nlev)           ! dry pressure level thickness     
     real(r8), intent(in)  :: ptop                               ! pressure at model top
     real(r8), intent(in)  :: p00                                ! reference pressure for Exner pressure (usually 1000hPa)
     real(r8), intent(in)  :: temp(i0:i1,j0:j1,nlev)             ! temperature
     real(r8), intent(out) :: theta_v(i0:i1,j0:j1,nlev)          ! virtual potential temperature

     real(r8) :: iexner(i0:i1,j0:j1,nlev)

     call get_exner(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                    dp_dry,ptop,p00,.true.,iexner)
     
     theta_v(:,:,:) = temp(:,:,:)*iexner(:,:,:)

   end subroutine get_virtual_theta
   !
   !****************************************************************************************************************   
   !   
   ! Compute geopotential from dry pressure level thichkness, water tracers, model top pressure and temperature
   !
   !****************************************************************************************************************   
   !   
   subroutine get_gz(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                     dp_dry,ptop,temp,phis,gz,pmid,dp,T_v)
     integer,  intent(in)  :: i0,i1,j0,j1,nlev,ntrac             ! array bounds
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)   ! tracer; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is dry mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(i0:i1,j0:j1,nlev)           ! dry pressure level thickness     
     real(r8), intent(in)  :: ptop                               ! pressure at model top
     real(r8), intent(in)  :: temp(i0:i1,j0:j1,nlev)             ! temperature
     real(r8), intent(in)  :: phis(i0:i1,j0:j1)                  ! surface geopotential
     real(r8), intent(out) :: gz(i0:i1,j0:j1,nlev)               ! geopotential
     real(r8), optional, intent(out) :: pmid(i0:i1,j0:j1,nlev)   ! mid-level pressure
     real(r8), optional, intent(out) :: dp(i0:i1,j0:j1,nlev)     ! pressure level thickness
     real(r8), optional, intent(out) :: t_v(i0:i1,j0:j1,nlev)    ! virtual temperature
     

     real(r8), dimension(i0:i1,j0:j1,nlev)   :: pmid_local, t_v_local, dp_local, R_dry
     real(r8), dimension(i0:i1,j0:j1,nlev+1) :: pint

     call get_pmid_from_dpdry(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                              dp_dry,ptop,pmid_local,pint=pint,dp=dp_local)       
     if (mixing_ratio==1) then 
       call get_virtual_temp(i0,i1,j0,j1,1,nlev,ntrac,tracer,t_v_local,temp=temp,&
            active_species_idx_dycore=active_species_idx)
       call get_R_dry(i0,i1,j0,j1,1,nlev,1,nlev,ntrac,tracer,active_species_idx,R_dry)
     else
       call get_virtual_temp(i0,i1,j0,j1,1,nlev,ntrac,tracer,t_v_local,temp=temp,dp_dry=dp_dry,&
            active_species_idx_dycore=active_species_idx)
       call get_R_dry(i0,i1,j0,j1,1,nlev,1,nlev,ntrac,tracer,active_species_idx, &
                      R_dry,fact=1.0_r8/dp_dry)
     end if
     call get_gz_given_dp_Tv_Rdry(i0,i1,j0,j1,nlev,dp_local,T_v_local,R_dry,phis,ptop,gz,pmid_local)

     if (present(pmid)) pmid=pmid_local
     if (present(T_v))  T_v=T_v_local
     if (present(dp))   dp=dp_local
   end subroutine get_gz
   !
   !****************************************************************************************************************   
   !   
   ! Compute geopotential from pressure level thickness and virtual temperature
   !
   !****************************************************************************************************************   
   !   
   subroutine get_gz_given_dp_Tv_Rdry(i0,i1,j0,j1,nlev,dp,T_v,R_dry,phis,ptop,gz,pmid)
     use dycore, only: dycore_is          
     integer,  intent(in)  :: i0,i1,j0,j1,nlev                 ! array bounds
     real(r8), intent(in)  :: dp   (i0:i1,j0:j1,nlev)          ! pressure level thickness 
     real(r8), intent(in)  :: T_v  (i0:i1,j0:j1,nlev)          ! virtual temperature
     real(r8), intent(in)  :: R_dry(i0:i1,j0:j1,nlev)          ! R dry
     real(r8), intent(in)  :: phis (i0:i1,j0:j1)               ! surface geopotential
     real(r8), intent(in)  :: ptop                             ! model top presure
     real(r8), intent(out) :: gz(i0:i1,j0:j1,nlev)             ! geopotential
     real(r8), optional, intent(out) :: pmid(i0:i1,j0:j1,nlev) ! mid-level pressure


     real(r8), dimension(i0:i1,j0:j1,nlev)   :: pmid_local
     real(r8), dimension(i0:i1,j0:j1,nlev+1) :: pint
     real(r8), dimension(i0:i1,j0:j1)        :: gzh, Rdry_tv
     integer :: k

     call get_pmid_from_dp(i0,i1,j0,j1,1,nlev,dp,ptop,pmid_local,pint)

     !
     ! integrate hydrostatic eqn
     !
     gzh = phis
     if (dycore_is ('LR')) then
       do k=nlev,1,-1
         Rdry_tv(:,:) = R_dry(:,:,k)*T_v(:,:,k)
         gz(:,:,k) = gzh(:,:)+Rdry_tv(:,:)*(1.0_r8-pint(:,:,k)/pmid_local(:,:,k))
         gzh(:,:)  = gzh(:,:) + Rdry_tv(:,:)*(log(pint(:,:,k+1))-log(pint(:,:,k)))
       end do
     else
       do k=nlev,1,-1
         Rdry_tv(:,:) = R_dry(:,:,k)*T_v(:,:,k)
         gz(:,:,k) = gzh(:,:)+Rdry_tv(:,:)*0.5_r8*dp(:,:,k)/pmid_local(:,:,k)
         gzh(:,:)  = gzh(:,:) + Rdry_tv(:,:)*dp(:,:,k)/pmid_local(:,:,k)
       end do       
     end if
     if (present(pmid)) pmid=pmid_local
   end subroutine get_gz_given_dp_Tv_Rdry
   !
   !****************************************************************************************************************   
   !   
   ! Compute Richardson number at cell interfaces (half levels)
   !
   !****************************************************************************************************************   
   !   
   subroutine get_Richardson_number(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx,&
        dp_dry,ptop,p00,temp,v,Richardson_number,pmid,dp)
     integer,  intent(in)  :: i0,i1,j0,j1,nlev,ntrac
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)     ! tracer; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                         ! 1 => tracer is dry mixing ratio
                                                                   ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)                ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(i0:i1,j0:j1,nlev)             ! dry pressure level thickness     
     real(r8), intent(in)  :: ptop                                 ! pressure at model top
     real(r8), intent(in)  :: p00                                  ! reference pressure for Exner pressure (usually 1000hPa)
     real(r8), intent(in)  :: temp(i0:i1,j0:j1,nlev)               ! temperature
     real(r8), intent(in)  :: v(i0:i1,j0:j1,2,nlev)                ! velocity components
     real(r8), intent(out) :: Richardson_number(i0:i1,j0:j1,nlev+1)!
     real(r8), optional, intent(out) :: pmid(i0:i1,j0:j1,nlev)   !
     real(r8), optional, intent(out) :: dp(i0:i1,j0:j1,nlev)   !     

     real(r8), dimension(i0:i1,j0:j1,nlev):: gz,theta_v
     real(r8), dimension(i0:i1,j0:j1)     :: pt1, pt2, phis
     integer :: k,km1
     real(r8), parameter:: ustar2 = 1.E-4_r8     

     phis = 0.0_r8
     if (present(pmid).and.present(dp)) then
        call get_gz(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                    dp_dry,ptop,temp,phis,gz,pmid=pmid,dp=dp)              
     else if (present(pmid)) then
        call get_gz(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                    dp_dry,ptop,temp,phis,gz,pmid=pmid)       
     else if (present(dp)) then
        call get_gz(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                    dp_dry,ptop,temp,phis,gz,dp=dp)       
     else
        call get_gz(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                    dp_dry,ptop,temp,phis,gz)
     end if
     call get_virtual_theta(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                            dp_dry,ptop,p00,temp,theta_v)
     Richardson_number(:,:,1)      = 0.0_r8
     Richardson_number(:,:,nlev+1) = 0.0_r8
     do k=nlev-1,2,-1
       km1=k-1
       pt1(:,:) = theta_v(:,:,km1)
       pt2(:,:) = theta_v(:,:,k)       
       Richardson_number(:,:,k) = (gz(:,:,km1)-gz(:,:,k))*(pt1-pt2)/( 0.5_r8*(pt1+pt2)*        &
            ((v(:,:,1,km1)-v(:,:,1,k))**2+(v(:,:,2,km1)-v(:,:,2,k))**2+ustar2) )
     end do
   end subroutine get_Richardson_number

   subroutine get_hydrostatic_static_energy(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx,&
        dp_dry,ptop,temp,phis,v,KE,thermalE,gz)
     integer,  intent(in)  :: i0,i1,j0,j1,nlev,ntrac
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)     ! tracer; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                         ! 1 => tracer is dry mixing ratio
                                                                   ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)                ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(i0:i1,j0:j1,nlev)             ! dry pressure level thickness     
     real(r8), intent(in)  :: ptop                                 ! pressure at model top
     real(r8), intent(in)  :: phis(i0:i1,j0:j1)                    ! surface geopotential
     real(r8), intent(in)  :: temp(i0:i1,j0:j1,nlev)               ! temperature
     real(r8), intent(in)  :: v(i0:i1,j0:j1,2,nlev)                ! velocity components

     real(r8), intent(out) :: KE(i0:i1,j0:j1,nlev),thermalE(i0:i1,j0:j1,nlev),gz(i0:i1,j0:j1,nlev)

     real(r8), dimension(i0:i1,j0:j1,nlev):: T_v,cp_dry

     call get_gz(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
                 dp_dry,ptop,temp,phis,gz,T_v=T_v)
     if (mixing_ratio==1) then
       call get_cp_dry(i0,i1,j0,j1,1,nlev,1,nlev,ntrac,tracer,active_species_idx,cp_dry)       
     else
        call get_cp_dry(i0,i1,j0,j1,1,nlev,1,nlev,ntrac,tracer,active_species_idx,cp_dry, &
                        fact=1.0_r8/dp_dry)              
     end if
     
     thermalE(:,:,:) = cp_dry(:,:,:)*T_v(:,:,:)
     KE(:,:,:)       = 0.5_r8*(v(:,:,2,:)**2+v(:,:,1,:)**2)
   end subroutine get_hydrostatic_static_energy
   !
   !****************************************************************************************************************   
   !
   ! get pressure from dry pressure and thermodynamic active species (e.g., forms of water: water vapor, cldliq, etc.)
   !
   !****************************************************************************************************************
   !     
   subroutine get_ps(i0,i1,j0,j1,k0,k1,ntrac,tracer_mass,active_species_idx,dp_dry,ps,ptop)
     integer,  intent(in)   :: i0,i1,j0,j1,k0,k1,ntrac
     real(r8), intent(in)   :: tracer_mass(i0:i1,j0:j1,k0:k1,1:ntrac) ! Tracer array
     real(r8), intent(in)   :: dp_dry(i0:i1,j0:j1,k0:k1)       ! dry pressure level thickness
     real(r8), intent(out)  :: ps(i0:i1,j0:j1)             ! surface pressure
     real(r8), intent(in)   :: ptop
     integer,  intent(in)   :: active_species_idx(:)

     integer                    :: i,j,k,m_cnst,nq
     real(r8)                   :: dp(i0:i1,j0:j1,k0:k1)       ! dry pressure level thickness
     
     dp = dp_dry
     do nq=dry_air_species_num+1,thermodynamic_active_species_num                   
       m_cnst = active_species_idx(nq)       
       do k=k0,k1
         do j=j0,j1
           do i = i0,i1
             dp(i,j,k) = dp(i,j,k) + tracer_mass(i,j,k,m_cnst)
           end do
         end do
       end do
     end do
     ps = ptop
     do k=k0,k1
       do j=j0,j1
         do i = i0,i1
           ps(i,j) = ps(i,j)+dp(i,j,k)
         end do
       end do
     end do
   end subroutine get_ps
   !
   !****************************************************************************************************************
   !  
   ! Compute dry air heaet capacity under constant pressure
   !
   !****************************************************************************************************************
   !     
   subroutine get_cp_dry(i0,i1,j0,j1,k0,k1,k0_trac,k1_trac,ntrac,tracer,active_species_idx,cp_dry,fact)
     integer,  intent(in)  :: i0,i1,j0,j1,k0,k1,ntrac,k0_trac,k1_trac
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,k0_trac:k1_trac,1:ntrac) ! Tracer array
     integer,  intent(in)  :: active_species_idx(:)
     real(r8), optional, intent(in) :: fact(i0:i1,j0:j1,k0_trac:k1_trac)       ! dry pressure level thickness
     real(r8), intent(out) :: cp_dry(i0:i1,j0:j1,k0:k1)       ! dry pressure level thickness     
     
     integer  :: i,j,k,m_cnst,nq
     real(r8) :: factor(i0:i1,j0:j1,k0_trac:k1_trac)       ! dry pressure level thickness          
     real(r8) :: residual(i0:i1,j0:j1,k0:k1), mm
     !
     ! dry air not species dependent
     !
     if (dry_air_species_num==0) then
       cp_dry = cpair
     else
       if (present(fact)) then
         factor = fact(:,:,:)
       else
         factor = 1.0_r8
       endif       

       cp_dry = 0.0_r8
       residual = 1.0_r8
       do nq=1,dry_air_species_num-1
         m_cnst = active_species_idx(nq)       
         do k=k0,k1
           do j=j0,j1
             do i = i0,i1
               mm = tracer(i,j,k,m_cnst)*factor(i,j,k)
               cp_dry(i,j,k) = cp_dry(i,j,k)+thermodynamic_active_species_cp(nq)*mm
               residual(i,j,k) = residual(i,j,k) - mm
             end do
           end do
         end do
       end do
       nq = dry_air_species_num
       do k=k0,k1
         do j=j0,j1
           do i = i0,i1
             cp_dry(i,j,k) = cp_dry(i,j,k)+thermodynamic_active_species_cp(nq)*residual(i,j,k)
           end do
         end do
       end do
     end if
   end subroutine get_cp_dry
   !
   !****************************************************************************************************************
   !  
   ! Compute generalized dry air gas constant R
   !
   !****************************************************************************************************************
   !     
   subroutine get_R_dry(i0,i1,j0,j1,k0,k1,k0_trac,k1_trac,ntrac,tracer,active_species_idx_dycore,R_dry,fact)
     integer,  intent(in)  :: i0,i1,j0,j1,k0,k1,ntrac,k0_trac,k1_trac    !array boundas
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,k0_trac:k1_trac,1:ntrac)!tracer array
     integer,  intent(in)  :: active_species_idx_dycore(:)               !index of active species in tracer
     real(r8), intent(out) :: R_dry(i0:i1,j0:j1,k0:k1)                   !dry air R
     real(r8), optional, intent(in) :: fact(i0:i1,j0:j1,k0_trac:k1_trac) !factor for converting tracer to dry mixing ratio
     
     integer :: i,j,k,m_cnst,nq
     real(r8):: factor(i0:i1,j0:j1,k0_trac:k1_trac), residual(i0:i1,j0:j1,k0:k1), mm
     if (dry_air_species_num==0) then
       !
       ! dry air not species dependent
       !
       R_dry = rair
     else
       if (present(fact)) then
         factor = fact(:,:,:)
       else
         factor = 1.0_r8
       endif       

       R_dry = 0.0_r8
       residual = 1.0_r8
       do nq=1,dry_air_species_num-1
         m_cnst = active_species_idx_dycore(nq)       
         do k=k0,k1
           do j=j0,j1
             do i = i0,i1
               mm = tracer(i,j,k,m_cnst)*factor(i,j,k)
               R_dry(i,j,k) = R_dry(i,j,k)+thermodynamic_active_species_R(nq)*mm
               residual(i,j,k) = residual(i,j,k) - mm
             end do
           end do
         end do
       end do
       !
       ! last dry air constituent derived from the others
       !
       nq = dry_air_species_num
       do k=k0,k1
         do j=j0,j1
           do i = i0,i1
             R_dry(i,j,k) = R_dry(i,j,k)+thermodynamic_active_species_R(nq)*residual(i,j,k)
           end do
         end do
       end do
     end if
   end subroutine get_R_dry
   !
   !*************************************************************************************************************************
   !
   ! Compute generalized R
   !
   !*************************************************************************************************************************
   !
   subroutine get_R(i0,i1,j0,j1,k0,k1,k0_trac,k1_trac,ntrac,tracer,active_species_idx,R,fact)
     integer,  intent(in)  :: i0,i1,j0,j1,k0,k1,ntrac,k0_trac,k1_trac
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,k0_trac:k1_trac,1:ntrac)!tracer array
     integer,  intent(in)  :: active_species_idx(:)                      !index of active species in tracer
     real(r8), intent(out) :: R(i0:i1,j0:j1,k0:k1)                       !generalized gas constant
     real(r8), optional, intent(in) :: fact(i0:i1,j0:j1,k0_trac:k1_trac) !factor for converting tracer to dry mixing ratio
     
     integer :: nq,itrac
     real(r8):: factor(i0:i1,j0:j1,k0_trac:k1_trac)
     real(r8), dimension(i0:i1,j0:j1,k0:k1)              :: sum_species
     integer, dimension(thermodynamic_active_species_num):: idx_local
     if (present(fact)) then
       call get_R_dry(i0,i1,j0,j1,k0,k1,k0_trac,k1_trac,ntrac,tracer,active_species_idx,R,fact=fact)
       factor = fact(:,:,:)
     else
       call get_R_dry(i0,i1,j0,j1,k0,k1,k0_trac,k1_trac,ntrac,tracer,active_species_idx,R)
       factor = 1.0_r8
     end if
     idx_local = active_species_idx
     sum_species = 1.0_r8 !all dry air species sum to 1
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = idx_local(nq)       
       sum_species(:,:,:) = sum_species(:,:,:) + tracer(:,:,:,itrac)*factor(:,:,:)
     end do
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = idx_local(nq)              
       R(:,:,:)      = R(:,:,:)+thermodynamic_active_species_R(nq)*tracer(:,:,:,itrac)*factor(:,:,:)
     end do
     R=R/sum_species
   end subroutine get_R
   !
   !*************************************************************************************************************************
   !
   ! compute molecular weight dry air
   !
   !*************************************************************************************************************************
   !
   subroutine get_mbarv(i0,i1,j0,j1,k0,k1,nlev,ntrac,tracer,active_species_idx,mbarv,fact)
     integer,  intent(in)  :: i0,i1,j0,j1,k0,k1,ntrac, nlev
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)   !tracer array
     integer,  intent(in)  :: active_species_idx(:)              !index of active species in tracer
     real(r8), intent(out) :: mbarv(i0:i1,j0:j1,k0:k1)           !molecular weight of dry air
     real(r8), optional, intent(in) :: fact(i0:i1,j0:j1,nlev)    !factor for converting tracer to dry mixing ratio
     
     integer :: i,j,k,m_cnst,nq
     real(r8):: factor(i0:i1,j0:j1,k0:k1), residual(i0:i1,j0:j1,k0:k1), mm
     !
     ! dry air not species dependent
     !
     if (dry_air_species_num==0) then
       mbarv = mwdry
     else
       if (present(fact)) then
         factor = fact(:,:,:)
       else
         factor = 1.0_r8
       endif       

       mbarv = 0.0_r8
       residual = 1.0_r8
       do nq=1,dry_air_species_num-1
         m_cnst = active_species_idx(nq)       
         do k=k0,k1
           do j=j0,j1
             do i = i0,i1
               mm = tracer(i,j,k,m_cnst)*factor(i,j,k)
               mbarv(i,j,k) = mbarv(i,j,k)+thermodynamic_active_species_mwi(nq)*mm
               residual(i,j,k) = residual(i,j,k) - mm
             end do
           end do
         end do
       end do
       nq = dry_air_species_num
       do k=k0,k1
         do j=j0,j1
           do i = i0,i1
             mbarv(i,j,k) = mbarv(i,j,k)+thermodynamic_active_species_mwi(nq)*residual(i,j,k)
           end do
         end do
       end do
       mbarv(i0:i1,j0:j1,k0:k1) = 1.0_r8/mbarv(i0:i1,j0:j1,k0:k1)
     end if
   end subroutine get_mbarv
   !
   !*************************************************************************************************************************
   !
   ! compute generalized kappa =Rdry/cpdry
   !
   !*************************************************************************************************************************
   !
   subroutine get_kappa_dry(i0,i1,j0,j1,k0,k1,nlev,ntrac,tracer,active_species_idx,kappa_dry,fact)
     integer,  intent(in)  :: i0,i1,j0,j1,k0,k1,ntrac,nlev
     real(r8), intent(in)  :: tracer(i0:i1,j0:j1,nlev,1:ntrac)   !tracer array
     integer,  intent(in)  :: active_species_idx(:)              !index of thermodynamic active tracers
     real(r8), intent(out) :: kappa_dry(i0:i1,j0:j1,k0:k1)       !kappa dry
     real(r8), optional, intent(in) :: fact(i0:i1,j0:j1,nlev)    !factor for converting tracer to dry mixing ratio
     !
     real(r8), allocatable, dimension(:,:,:) :: cp_dry,R_dry
     !
     ! dry air not species dependent
     if (dry_air_species_num==0) then
       kappa_dry= rair/cpair
     else
       allocate(R_dry(i0:i1,j0:j1,k0:k1))
       allocate(cp_dry(i0:i1,j0:j1,k0:k1))
       if (present(fact)) then
         call get_cp_dry(i0,i1,j0,j1,k0,k1,1,nlev,ntrac,tracer,active_species_idx,cp_dry,fact=fact)         
         call get_R_dry(i0,i1,j0,j1,k0,k1,1,nlev,ntrac,tracer,active_species_idx,R_dry,fact=fact)         
       else
         call get_cp_dry(i0,i1,j0,j1,k0,k1,1,nlev,ntrac,tracer,active_species_idx,cp_dry)         
         call get_R_dry(i0,i1,j0,j1,k0,k1,1,nlev,ntrac,tracer,active_species_idx,R_dry)         
       end if
       kappa_dry = R_dry/cp_dry       
       deallocate(R_dry,cp_dry)
     end if
   end subroutine get_kappa_dry
   !
   !****************************************************************************************************************
   !
   ! Compute sum of thermodynamically active species
   !
   ! tracer is in units of dry mixing ratio unless optional argument dp_dry is present in which case tracer is
   ! in units of "mass" (=m*dp)
   !
   !****************************************************************************************************************
   !
   subroutine get_sum_species(i0,i1,j0,j1,k0,k1,ntrac,tracer,active_species_idx,sum_species,dp_dry)
     integer,  intent(in)           :: i0,i1,j0,j1,k0,k1,ntrac
     real(r8), intent(in)           :: tracer(i0:i1,j0:j1,k0:k1,1:ntrac)   ! tracer array
     integer,  intent(in)           :: active_species_idx(:)               ! index for thermodynamic active tracers
     real(r8), optional, intent(in) :: dp_dry(i0:i1,j0:j1,k0:k1)           ! dry pressure level thickness is present
                                                                           ! then tracer is in units of mass
     real(r8), intent(out)          :: sum_species(i0:i1,j0:j1,k0:k1)      ! sum species

     real(r8) :: factor(i0:i1,j0:j1,k0:k1)
     integer  :: nq,itrac
     
     if (present(dp_dry)) then
       factor = 1.0_r8/dp_dry(:,:,:)
     else
       factor = 1.0_r8
     endif
     sum_species = 1.0_r8 !all dry air species sum to 1
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = active_species_idx(nq)
       sum_species(:,:,:) = sum_species(:,:,:) + tracer(:,:,:,itrac)*factor(:,:,:)
     end do
   end subroutine get_sum_species
   !
   !****************************************************************************************************************
   !
   ! g*compute thermal energy = cp*T*dp, where dp is pressure level thickness, cp is generalized cp and T temperature
   !
   ! Note:tracer is in units of m*dp_dry ("mass")
   !
   !****************************************************************************************************************   
   !
   subroutine get_thermal_energy(i0,i1,j0,j1,k0,k1,ntrac,tracer_mass,temp,dp_dry,thermal_energy, &
                                 active_species_idx_dycore)
     integer,  intent(in)           :: i0,i1,j0,j1,k0,k1,ntrac      
     real(r8), intent(in)           :: tracer_mass(i0:i1,j0:j1,k0:k1,ntrac)!tracer array (mass weighted)
     real(r8), intent(in)           :: temp(i0:i1,j0:j1,k0:k1)             !temperature
     real(r8), intent(in)           :: dp_dry(i0:i1,j0:j1,k0:k1)           !dry presure level thickness
     real(r8), optional, intent(out):: thermal_energy(i0:i1,j0:j1,k0:k1)   !thermal energy in each column: sum cp*T*dp
     !
     ! array of indicies for index of thermodynamic active species in dycore tracer array
     ! (if different from physics index)
     !
     integer, optional, dimension(:) :: active_species_idx_dycore
     
     ! local vars
     integer :: nq, itrac
     integer, dimension(thermodynamic_active_species_num)                   :: idx_local
     !
     ! some sanity checks
     !
     if (present(active_species_idx_dycore)) then
       idx_local = active_species_idx_dycore
     else
       idx_local = thermodynamic_active_species_idx
     end if     
     !
     ! "mass-weighted" cp (dp must be dry)
     !
     if (dry_air_species_num==0) then
       thermal_energy(:,:,:) = thermodynamic_active_species_cp(0)*dp_dry(:,:,:)
     else
       call get_cp_dry(i0,i1,j0,j1,k0,k1,k0,k1,ntrac,tracer_mass,idx_local,thermal_energy,fact=1.0_r8/dp_dry(:,:,:))
       thermal_energy(:,:,:) = thermal_energy(:,:,:)*dp_dry(:,:,:)
     end if
     !
     ! tracer is in units of m*dp ("mass"), where m is dry mixing ratio and dry pressure level thickness
     !         
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = idx_local(nq)
       thermal_energy(:,:,:) = thermal_energy(:,:,:)+thermodynamic_active_species_cp(nq)*tracer_mass(:,:,:,itrac)
     end do
     thermal_energy(:,:,:) = thermal_energy(:,:,:)*temp(:,:,:)
   end subroutine get_thermal_energy
   !
   !****************************************************************************************************************
   !  
   ! Compute virtual temperature T_v
   !
   ! tracer is in units of dry mixing ratio unless optional argument dp_dry is present in which case tracer is
   ! in units of "mass" (=m*dp)
   !
   ! If temperature is not supplied then just return factor that T needs to be multiplied by to get T_v
   !
   !****************************************************************************************************************
   !
   subroutine get_virtual_temp(i0,i1,j0,j1,k0,k1,ntrac,tracer,T_v,temp,dp_dry,sum_q, &
                               active_species_idx_dycore)
     use cam_logfile,     only: iulog
     ! args
     integer,  intent(in)           :: i0,i1,j0,j1,k0,k1,ntrac
     real(r8), intent(in)           :: tracer(i0:i1,j0:j1,k0:k1,ntrac) !tracer array
     real(r8), intent(out)          :: T_v(i0:i1,j0:j1,k0:k1)          !virtual temperature
     real(r8), optional, intent(in) :: temp(i0:i1,j0:j1,k0:k1)         !temperature
     real(r8), optional, intent(in) :: dp_dry(i0:i1,j0:j1,k0:k1)       !dry pressure level thickness
     real(r8), optional,intent(out) :: sum_q(i0:i1,j0:j1,k0:k1)        !sum tracer
     !
     ! array of indicies for index of thermodynamic active species in dycore tracer array
     ! (if different from physics index)
     !
     integer, optional,  intent(in) :: active_species_idx_dycore(:)
     
     ! local vars
     integer :: itrac,nq
     real(r8),  dimension(i0:i1,j0:j1,k0:k1)              :: sum_species, factor, Rd
     integer, dimension(thermodynamic_active_species_num) :: idx_local,idx

     if (present(active_species_idx_dycore)) then
       idx_local = active_species_idx_dycore
     else
       idx_local = thermodynamic_active_species_idx
     end if

     if (present(dp_dry)) then
       factor = 1.0_r8/dp_dry
     else
       factor = 1.0_r8       
     end if

     sum_species = 1.0_r8 !all dry air species sum to 1
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = idx_local(nq)       
       sum_species(:,:,:) = sum_species(:,:,:) + tracer(:,:,:,itrac)*factor(:,:,:)
     end do
   
     call get_R_dry (i0,i1,j0,j1,k0,k1,k0,k1,ntrac,tracer,idx_local,Rd,fact=factor)
     t_v(:,:,:)  = Rd(:,:,:)
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = idx_local(nq)
       t_v(:,:,:) = t_v(:,:,:)+thermodynamic_active_species_R(nq)*tracer(:,:,:,itrac)*factor(:,:,:)
     end do
     if (present(temp)) then
       t_v(:,:,:)  = t_v(:,:,:)*temp(:,:,:)/(Rd(:,:,:)*sum_species)
     else
       t_v(:,:,:)  = t_v(:,:,:)/(Rd(:,:,:)*sum_species)       
     end if
     if (present(sum_q)) sum_q=sum_species
   end subroutine get_virtual_temp
   !
   !*************************************************************************************************************************
   !
   ! Compute generalized heat capacity at constant pressure
   !
   !*************************************************************************************************************************
   !
   subroutine get_cp(i0,i1,j0,j1,k0,k1,ntrac,tracer,inv_cp,cp,dp_dry,active_species_idx_dycore)
     use cam_logfile,     only: iulog
     ! args
     integer,  intent(in)           :: i0,i1,j0,j1,k0,k1,ntrac
     real(r8), intent(in)           :: tracer(i0:i1,j0:j1,k0:k1,ntrac) ! Tracer array
     real(r8), optional, intent(in) :: dp_dry(i0:i1,j0:j1,k0:k1)
     logical , intent(in)           :: inv_cp!output inverse cp instead of cp
     real(r8), intent(out)          :: cp(i0:i1,j0:j1,k0:k1)
     !
     ! array of indicies for index of thermodynamic active species in dycore tracer array
     ! (if different from physics index)
     !
     integer, optional, intent(in)  :: active_species_idx_dycore(:)
     
     ! local vars
     integer :: nq,i,j,k, itrac
     real(r8),  dimension(i0:i1,j0:j1,k0:k1)  :: sum_species, sum_cp, factor
     integer, dimension(thermodynamic_active_species_num) :: idx_local

     if (present(active_species_idx_dycore)) then
       idx_local = active_species_idx_dycore
     else
       idx_local = thermodynamic_active_species_idx
     end if

     if (present(dp_dry)) then
       factor = 1.0_r8/dp_dry
     else
       factor = 1.0_r8       
     end if

     sum_species = 1.0_r8 !all dry air species sum to 1
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = idx_local(nq)       
       sum_species(:,:,:) = sum_species(:,:,:) + tracer(:,:,:,itrac)*factor(:,:,:)
     end do

     if (dry_air_species_num==0) then
       sum_cp = thermodynamic_active_species_cp(0)
     else
       call get_cp_dry(i0,i1,j0,j1,k0,k1,k0,k1,ntrac,tracer,idx_local,sum_cp,fact=factor)
     end if
     do nq=dry_air_species_num+1,thermodynamic_active_species_num
       itrac = idx_local(nq)              
       sum_cp(:,:,:)      = sum_cp(:,:,:)+thermodynamic_active_species_cp(nq)*tracer(:,:,:,itrac)*factor(:,:,:)
     end do
     if (inv_cp) then
       cp=sum_species/sum_cp
     else
       cp=sum_cp/sum_species
     end if

   end subroutine get_cp
   !
   !*************************************************************************************************************************
   !
   ! compute reference pressure levels
   !
   !*************************************************************************************************************************
   !
   subroutine get_dp_ref(hyai, hybi, ps0, i0,i1,j0,j1,k0,k1,phis,dp_ref,ps_ref)
     integer,  intent(in)           :: i0,i1,j0,j1,k0,k1
     real(r8), intent(in)           :: hyai(k0:k1+1),hybi(k0:k1+1),ps0
     real(r8), intent(in)           :: phis(i0:i1,j0:j1)
     real(r8), intent(out)          :: dp_ref(i0:i1,j0:j1,k0:k1)
     real(r8), intent(out)          :: ps_ref(i0:i1,j0:j1)
     integer :: k
     !
     ! use static reference pressure (hydrostatic balance incl. effect of topography)
     !
     ps_ref(:,:) = ps0*exp(-phis(:,:)/(Rair*Tref))
     do k=k0,k1
       dp_ref(:,:,k) = ((hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_ref(:,:))
     end do
   end subroutine get_dp_ref
   !
   !*************************************************************************************************************************
   !
   ! compute dry densisty from temperature (temp) and pressure (dp_dry and tracer)
   !
   !*************************************************************************************************************************
   !
   subroutine get_rho_dry(i0,i1,j0,j1,k1,nlev,ntrac,tracer,temp,ptop,dp_dry,tracer_mass,&
        rho_dry, rhoi_dry,active_species_idx_dycore,pint_out,pmid_out)
     ! args
     integer,  intent(in)           :: i0,i1,j0,j1,k1,ntrac,nlev
     real(r8), intent(in)           :: tracer(i0:i1,j0:j1,nlev,ntrac) ! Tracer array
     real(r8), intent(in)           :: temp(i0:i1,j0:j1,1:nlev) ! Temperature
     real(r8), intent(in)           :: ptop
     real(r8), intent(in)           :: dp_dry(i0:i1,j0:j1,nlev)
     logical,  intent(in)           :: tracer_mass
     real(r8), optional,intent(out) :: rho_dry(i0:i1,j0:j1,1:k1)
     real(r8), optional,intent(out) :: rhoi_dry(i0:i1,j0:j1,1:k1+1)
     !
     ! array of indicies for index of thermodynamic active species in dycore tracer array
     ! (if different from physics index)
     !
     integer, optional, intent(in)   :: active_species_idx_dycore(:)
     real(r8),optional,intent(out)   :: pint_out(i0:i1,j0:j1,1:k1+1)
     real(r8),optional,intent(out)   :: pmid_out(i0:i1,j0:j1,1:k1)
     
     ! local vars
     integer :: i,j,k
     real(r8),  dimension(i0:i1,j0:j1,1:k1)              :: pmid
     real(r8):: pint(i0:i1,j0:j1,1:k1+1)
     real(r8), allocatable :: R_dry(:,:,:)
     integer,  dimension(thermodynamic_active_species_num):: idx_local

     if (present(active_species_idx_dycore)) then
       idx_local = active_species_idx_dycore
     else
       idx_local = thermodynamic_active_species_idx
     end if
     !
     ! we assume that air is dry where molecular viscosity may be significant
     !
     call get_pmid_from_dp(i0,i1,j0,j1,1,k1,dp_dry,ptop,pmid,pint=pint)
     if (present(pint_out)) pint_out=pint
     if (present(pint_out)) pmid_out=pmid
     if (present(rhoi_dry)) then
       allocate(R_dry(i0:i1,j0:j1,1:k1+1))
       if (tracer_mass) then           
         call get_R_dry(i0,i1,j0,j1,1,k1+1,1,nlev,ntrac,tracer,idx_local,R_dry,fact=1.0_r8/dp_dry)
       else
         call get_R_dry(i0,i1,j0,j1,1,k1+1,1,nlev,ntrac,tracer,idx_local,R_dry)
       end if
       do k=2,k1+1
         rhoi_dry(i0:i1,j0:j1,k) = 0.5_r8*(temp(i0:i1,j0:j1,k)+temp(i0:i1,j0:j1,k-1))!could be more accurate!
         rhoi_dry(i0:i1,j0:j1,k) = pint(i0:i1,j0:j1,k)/(rhoi_dry(i0:i1,j0:j1,k)*R_dry(i0:i1,j0:j1,k)) !ideal gas law for dry air
       end do
       !
       ! extrapolate top level value
       !
       k=1
       rhoi_dry(i0:i1,j0:j1,k) = 1.5_r8*(temp(i0:i1,j0:j1,1)-0.5_r8*temp(i0:i1,j0:j1,2))
       rhoi_dry(i0:i1,j0:j1,k) = pint(i0:i1,j0:j1,1)/(rhoi_dry(i0:i1,j0:j1,k)*R_dry(i0:i1,j0:j1,k)) !ideal gas law for dry air
       deallocate(R_dry)
     end if
     if (present(rho_dry)) then
       allocate(R_dry(i0:i1,j0:j1,1:k1))
       if (tracer_mass) then
         call get_R_dry(i0,i1,j0,j1,1,k1,1,nlev,ntrac,tracer,idx_local,R_dry,fact=1.0_r8/dp_dry)
       else
         call get_R_dry(i0,i1,j0,j1,1,k1,1,nlev,ntrac,tracer,idx_local,R_dry)
       end if
       do k=1,k1
         do j=j0,j1
           do i=i0,i1
             rho_dry(i,j,k) = pmid(i,j,k)/(temp(i,j,k)*R_dry(i,j,k)) !ideal gas law for dry air
           end do
         end do
       end do
     end if
   end subroutine get_rho_dry
   !
   !*************************************************************************************************************************
   !
   ! compute 3D molecular diffusion and thermal conductivity
   !   
   !*************************************************************************************************************************
   !
   subroutine get_molecular_diff_coef(i0,i1,j0,j1,k1,nlev,temp,get_at_interfaces,sponge_factor,kmvis,kmcnd, ntrac,&
        tracer, fact, active_species_idx_dycore, mbarv_in)
     ! args
     integer,  intent(in)           :: i0,i1,j0,j1,k1,nlev
     real(r8), intent(in)           :: temp(i0:i1,j0:j1,nlev) ! temperature
     integer,  intent(in)           :: get_at_interfaces      ! 1: compute kmvis and kmcnd at interfaces
                                                              ! 0: compute kmvis and kmcnd at mid-levels
     real(r8), intent(in)           :: sponge_factor(1:k1)    ! multiply kmvis and kmcnd with sponge_factor (for sponge layer)
     real(r8), intent(out)          :: kmvis(i0:i1,j0:j1,1:k1+get_at_interfaces)
     real(r8), intent(out)          :: kmcnd(i0:i1,j0:j1,1:k1+get_at_interfaces)
     integer , intent(in)           :: ntrac
     real(r8), intent(in)           :: tracer(i0:i1,j0:j1,nlev,1:ntrac)   ! tracer array
     integer,  intent(in), optional :: active_species_idx_dycore(:)  ! index of active species in tracer
     real(r8), intent(in), optional :: fact(i0:i1,j0:j1,k1)          ! if tracer is in units of mass or moist
                                                                          ! fact converts to dry mixing ratio: tracer/fact
     real(r8), intent(in), optional :: mbarv_in(i0:i1,j0:j1,1:k1)         ! composition dependent atmosphere mean mass
     !
     ! local vars
     !
     integer :: i,j,k,icnst,ispecies
     real(r8):: mbarvi,mm,residual             ! Mean mass at mid level
     real(r8):: cnst_vis, cnst_cnd, temp_local
     real(r8), dimension(i0:i1,j0:j1,1:k1)                :: factor,mbarv       
     integer,  dimension(thermodynamic_active_species_num):: idx_local

     !--------------------------------------------
     ! Set constants needed for updates
     !--------------------------------------------

     if (dry_air_species_num==0) then

       cnst_vis = (kv1*mmro2*o2_mwi + kv2*mmrn2*n2_mwi)*mbar*1.e-7_r8
       cnst_cnd = (kc1*mmro2*o2_mwi + kc2*mmrn2*n2_mwi)*mbar*1.e-5_r8
       if (get_at_interfaces==1) then
           do k=2,k1
             do j=j0,j1
               do i=i0,i1
                 temp_local   = 0.5_r8*(temp(i,j,k)+temp(i,j,k-1))
                 kmvis(i,j,k) = sponge_factor(k)*cnst_vis*temp_local**kv4                               
                 kmcnd(i,j,k) = sponge_factor(k)*cnst_cnd*temp_local**kc4
               end do
             end do
           end do
           !
           ! extrapolate top level value
           !
           kmvis(i0:i1,j0:j1,1) = 1.5_r8*kmvis(i0:i1,j0:j1,2)-0.5_r8*kmvis(i0:i1,j0:j1,3)
           kmcnd(i0:i1,j0:j1,1) = 1.5_r8*kmcnd(i0:i1,j0:j1,2)-0.5_r8*kmcnd(i0:i1,j0:j1,3)
       else if (get_at_interfaces==0) then 
         do k=1,k1
           do j=j0,j1
             do i=i0,i1
               kmvis(i,j,k) = sponge_factor(k)*cnst_vis*temp(i,j,k)**kv4                               
               kmcnd(i,j,k) = sponge_factor(k)*cnst_cnd*temp(i,j,k)**kc4
             end do
           end do
         end do
       else
         call endrun('get_molecular_diff_coef: get_at_interfaces must be 0 or 1')
       end if
     else
       if (present(active_species_idx_dycore)) then
         idx_local = active_species_idx_dycore
       else
         idx_local = thermodynamic_active_species_idx
       end if
       if (present(fact)) then
         factor = fact(:,:,:)
       else
         factor = 1.0_r8
       endif
       if (present(mbarv_in)) then
         mbarv = mbarv_in
       else
         call get_mbarv(i0,i1,j0,j1,1,k1,nlev,ntrac,tracer,idx_local,mbarv,fact=factor)
       end if
       !
       ! major species dependent code
       !
       if (get_at_interfaces==1) then
         do k=2,k1
           do j=j0,j1
             do i=i0,i1               
               kmvis(i,j,k) = 0.0_r8
               kmcnd(i,j,k) = 0.0_r8
               residual = 1.0_r8
               do icnst=1,dry_air_species_num-1             
                 ispecies = idx_local(icnst)
                 mm       = 0.5_r8*(tracer(i,j,k,ispecies)*factor(i,j,k)+tracer(i,j,k-1,ispecies)*factor(i,j,k-1))
                 kmvis(i,j,k) = kmvis(i,j,k)+thermodynamic_active_species_kv(icnst)* &
                                             thermodynamic_active_species_mwi(icnst)*mm
                 kmcnd(i,j,k) = kmcnd(i,j,k)+thermodynamic_active_species_kc(icnst)* &
                                             thermodynamic_active_species_mwi(icnst)*mm
                 residual         = residual - mm
               end do
               icnst=dry_air_species_num
               ispecies = idx_local(icnst)
               kmvis(i,j,k) = kmvis(i,j,k)+thermodynamic_active_species_kv(icnst)* &
                                           thermodynamic_active_species_mwi(icnst)*residual
               kmcnd(i,j,k) = kmcnd(i,j,k)+thermodynamic_active_species_kc(icnst)* &
                                           thermodynamic_active_species_mwi(icnst)*residual
               
               temp_local = .5_r8*(temp(i,j,k-1)+temp(i,j,k))
               mbarvi = 0.5_r8*(mbarv(i,j,k-1)+mbarv(i,j,k))
               kmvis(i,j,k) = kmvis(i,j,k)*mbarvi*temp_local**kv4*1.e-7_r8                              
               kmcnd(i,j,k) = kmcnd(i,j,k)*mbarvi*temp_local**kc4*1.e-5_r8
             enddo
           enddo
         end do
         do j=j0,j1
           do i=i0,i1
             kmvis(i,j,1)    = 1.5_r8*kmvis(i,j,2)-.5_r8*kmvis(i,j,3)
             kmcnd(i,j,1)    = 1.5_r8*kmcnd(i,j,2)-.5_r8*kmcnd(i,j,3)
             kmvis(i,j,k1+1) = kmvis(i,j,k1)
             kmcnd(i,j,k1+1) = kmcnd(i,j,k1)
           end do
         end do
       else if (get_at_interfaces==0) then 
       else
         call endrun('get_molecular_diff_coef: get_at_interfaces must be 0 or 1')
       end if
     end if
   end subroutine get_molecular_diff_coef
   !
   !*************************************************************************************************************************
   !
   ! compute reference vertical profile of density, molecular diffusion and thermal conductivity
   !   
   !*************************************************************************************************************************
   !
   subroutine get_molecular_diff_coef_reference(k0,k1,tref,press,sponge_factor,kmvis_ref,kmcnd_ref,rho_ref)
     ! args
     integer,  intent(in)           :: k0,k1                !min/max vertical index
     real(r8), intent(in)           :: tref                 !reference temperature
     real(r8), intent(in)           :: press(k0:k1)         !pressure
     real(r8), intent(in)           :: sponge_factor(k0:k1) !multiply kmvis and kmcnd with sponge_factor (for sponge layer)
     real(r8), intent(out)          :: kmvis_ref(k0:k1)     !reference molecular diffusion coefficient
     real(r8), intent(out)          :: kmcnd_ref(k0:k1)     !reference thermal conductivity coefficient
     real(r8), intent(out)          :: rho_ref(k0:k1)       !reference density
     
     ! local vars
     integer :: k
    
     !--------------------------------------------
     ! Set constants needed for updates
     !--------------------------------------------

     do k=k0,k1
       rho_ref(k) = press(k)/(tref*Rair) !ideal gas law for dry air
       kmvis_ref(k) = sponge_factor(k)* &
            (kv1*mmro2*o2_mwi +         &
             kv2*mmrn2*n2_mwi)*mbar*    &
             tref**kv4 * 1.e-7_r8
       kmcnd_ref(k) = sponge_factor(k)* &
            (kc1*mmro2*o2_mwi +         &
             kc2*mmrn2*n2_mwi)*mbar*    &
             tref**kc4 * 1.e-5_r8
     end do
   end subroutine get_molecular_diff_coef_reference


end module physconst
