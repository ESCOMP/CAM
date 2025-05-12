! cam_thermo module provides interfaces to compute thermodynamic quantities
module cam_thermo

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use cam_abortutils,  only: endrun
   use air_composition, only: thermodynamic_active_species_num
   use air_composition, only: thermodynamic_active_species_idx
   use air_composition, only: thermodynamic_active_species_idx_dycore
   use air_composition, only: thermodynamic_active_species_cp
   use air_composition, only: thermodynamic_active_species_R
   use air_composition, only: thermodynamic_active_species_mwi
   use air_composition, only: thermodynamic_active_species_kv
   use air_composition, only: thermodynamic_active_species_kc
   use air_composition, only: thermodynamic_active_species_liq_num
   use air_composition, only: thermodynamic_active_species_ice_num
   use air_composition, only: thermodynamic_active_species_liq_idx
   use air_composition, only: thermodynamic_active_species_liq_idx_dycore
   use air_composition, only: thermodynamic_active_species_ice_idx
   use air_composition, only: thermodynamic_active_species_ice_idx_dycore
   use air_composition, only: dry_air_species_num
   use air_composition, only: enthalpy_reference_state
   use air_composition, only: mmro2, mmrn2, o2_mwi, n2_mwi, mbar

   implicit none
   private
   save

   ! subroutines to compute thermodynamic quantities
   !
   ! See Lauritzen et al. (2018) for formulae
   !     DOI: 10.1029/2017MS001257
   !     https://opensky.ucar.edu/islandora/object/articles:21929

   ! cam_thermo_init: Initialize constituent dependent properties
   public :: cam_thermo_init
   ! cam_thermo_dry_air_update: Update dry air composition dependent properties
   public :: cam_thermo_dry_air_update
   ! cam_thermo_water_update: Update water dependent properties
   public :: cam_thermo_water_update
   ! get_enthalpy: enthalpy quantity = dp*cp*T
   public :: get_enthalpy
   ! get_virtual_temp: virtual temperature
   public :: get_virtual_temp
   ! get_sum_species: sum of thermodynamically active species:
   !                  Note: dp = dp_dry * sum_species
   public :: get_sum_species
   ! get_virtual_theta: virtual potential temperature
   public :: get_virtual_theta
   ! cam_thermo_calc_kappav: update species dependent kappa for FV dycore
   public :: cam_thermo_calc_kappav
   ! get_dp: pressure level thickness from dry dp and dry mixing ratios
   public :: get_dp
   ! get_pmid_from_dp: full level pressure from dp (approximation depends on dycore)
   public :: get_pmid_from_dp
   ! get_ps: surface pressure
   public :: get_ps
   ! get_gz: geopotential
   public :: get_gz
   ! get_Richardson_number: Richardson number at layer interfaces
   public :: get_Richardson_number
   ! get_kappa_dry: (generalized) dry kappa = R_dry/cp_dry
   public :: get_kappa_dry
   ! get_dp_ref: reference pressure layer thickness (include topography)
   public :: get_dp_ref
   ! get_molecular_diff_coef: molecular diffusion and thermal conductivity
   public :: get_molecular_diff_coef
   ! get_molecular_diff_coef_reference: reference vertical profile of density,
   !                             molecular diffusion and thermal conductivity
   public :: get_molecular_diff_coef_reference
   ! get_rho_dry: dry density from temperature (temp) and
   !              pressure (dp_dry and tracer)
   public :: get_rho_dry
   ! get_exner: Exner pressure
   public :: get_exner
   ! get_hydrostatic_energy: Vertically integrated total energy
   public :: get_hydrostatic_energy

   ! Public variables
   ! mixing_ratio options
   integer, public, parameter :: DRY_MIXING_RATIO = 1
   integer, public, parameter :: MASS_MIXING_RATIO = 2
   !---------------  Variables below here are for WACCM-X ---------------------
   ! kmvis: molecular viscosity      kg/m/s
   real(r8), public, protected, allocatable :: kmvis(:,:,:)
   ! kmcnd: molecular conductivity   J/m/s/K
   real(r8), public, protected, allocatable :: kmcnd(:,:,:)

   !------------- Variables for consistent themodynamics --------------------
   !

   !
   ! Interfaces for public routines
   interface get_gz
      ! get_gz_geopotential (with dp_dry, ptop, temp, and phis as input)
      module procedure get_gz_from_dp_dry_ptop_temp_1hd
      ! get_gz_given_dp_Tv_Rdry: geopotential (with dp,dry R and Tv as input)
      module procedure get_gz_given_dp_Tv_Rdry_1hd
      module procedure get_gz_given_dp_Tv_Rdry_2hd
   end interface get_gz

   interface get_enthalpy
      module procedure get_enthalpy_1hd
      module procedure get_enthalpy_2hd
   end interface get_enthalpy

   interface get_virtual_temp
      module procedure get_virtual_temp_1hd
      module procedure get_virtual_temp_2hd
   end interface get_virtual_temp

   interface get_sum_species
      module procedure get_sum_species_1hd
      module procedure get_sum_species_2hd
   end interface get_sum_species

   interface get_dp
      module procedure get_dp_1hd
      module procedure get_dp_2hd
   end interface get_dp

   interface get_pmid_from_dp
      module procedure get_pmid_from_dpdry_1hd
      module procedure get_pmid_from_dp_1hd
   end interface get_pmid_from_dp

   interface get_exner
      module procedure get_exner_1hd
   end interface get_exner

   interface get_virtual_theta
      module procedure get_virtual_theta_1hd
   end interface get_virtual_theta

   interface get_Richardson_number
      module procedure get_Richardson_number_1hd
   end interface get_Richardson_number

   interface get_ps
      module procedure get_ps_1hd
      module procedure get_ps_2hd
   end interface get_ps

   interface get_kappa_dry
      module procedure get_kappa_dry_1hd
      module procedure get_kappa_dry_2hd
   end interface get_kappa_dry

   interface get_dp_ref
      module procedure get_dp_ref_1hd
      module procedure get_dp_ref_2hd
   end interface get_dp_ref

   interface get_rho_dry
      module procedure get_rho_dry_1hd
      module procedure get_rho_dry_2hd
   end interface get_rho_dry

   interface get_molecular_diff_coef
      module procedure get_molecular_diff_coef_1hd
      module procedure get_molecular_diff_coef_2hd
   end interface get_molecular_diff_coef

   interface cam_thermo_calc_kappav
      ! Since this routine is currently only used by the FV dycore,
      !    a 1-d interface is not needed (but can easily be added)
      module procedure cam_thermo_calc_kappav_2hd
   end interface cam_thermo_calc_kappav

   interface get_hydrostatic_energy
      module procedure get_hydrostatic_energy_1hd
      ! This routine is currently only called from the physics so a
      !    2-d interface is not needed (but can easily be added)
   end interface get_hydrostatic_energy

   integer, public, parameter :: thermo_budget_num_vars = 10
   integer, public, parameter :: wvidx = 1
   integer, public, parameter :: wlidx = 2
   integer, public, parameter :: wiidx = 3
   integer, public, parameter :: seidx = 4 ! enthalpy or internal energy (W/m2) index
   integer, public, parameter :: poidx = 5 ! surface potential or potential energy index
   integer, public, parameter :: keidx = 6 ! kinetic energy index
   integer, public, parameter :: mridx = 7
   integer, public, parameter :: moidx = 8
   integer, public, parameter :: ttidx = 9
   integer, public, parameter :: teidx = 10
   character (len = 2)  ,public, dimension(thermo_budget_num_vars) :: thermo_budget_vars  = &
        (/"WV"  ,"WL"  ,"WI"  ,"SE"   ,"PO"   ,"KE"   ,"MR"   ,"MO"   ,"TT"   ,"TE"   /)
   character (len = 46) ,public, dimension(thermo_budget_num_vars) :: thermo_budget_vars_descriptor = (/&
        "Total column water vapor                      ",&
        "Total column liquid water                     ",&
        "Total column frozen water                     ",&
        "Total column enthalpy or internal energy      ",&
        "Total column srf potential or potential energy",&
        "Total column kinetic energy                   ",&
        "Total column wind axial angular momentum      ",&
        "Total column mass axial angular momentum      ",&
        "Total column test_tracer                      ",&
        "Total column energy (ke + se + po)            "/)

   character (len = 14), public, dimension(thermo_budget_num_vars)  :: &
        thermo_budget_vars_unit = (/&
        "kg/m2        ","kg/m2        ","kg/m2        ","J/m2         ",&
        "J/m2         ","J/m2         ","kg*m2/s*rad2 ","kg*m2/s*rad2 ",&
        "kg/m2        ","J/m2         "/)
   logical              ,public, dimension(thermo_budget_num_vars) :: thermo_budget_vars_massv = (/&
        .true.,.true.,.true.,.false.,.false.,.false.,.false.,.false.,.true.,.false./)
CONTAINS

   !===========================================================================

   subroutine cam_thermo_init()
      use shr_infnan_mod,  only: assignment(=), shr_infnan_qnan
      use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk

      integer                     :: ierr
      character(len=*), parameter :: subname = "cam_thermo_init"
      character(len=*), parameter :: errstr = subname//": failed to allocate "

      !------------------------------------------------------------------------
      !  Allocate constituent dependent properties
      !------------------------------------------------------------------------
      allocate(kmvis(pcols,pverp,begchunk:endchunk), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"kmvis")
      end if
      allocate(kmcnd(pcols,pverp,begchunk:endchunk), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"kmcnd")
      end if

      !------------------------------------------------------------------------
      !  Initialize constituent dependent properties
      !------------------------------------------------------------------------
      kmvis(:pcols,  :pver, begchunk:endchunk) = shr_infnan_qnan
      kmcnd(:pcols,  :pver, begchunk:endchunk) = shr_infnan_qnan

   end subroutine cam_thermo_init
   !
   !***************************************************************************
   !
   ! cam_thermo_dry_air_update: update dry air species dependent constants for physics
   !
   !***************************************************************************
   !
   subroutine cam_thermo_dry_air_update(mmr, T, lchnk, ncol, to_dry_factor)
      use air_composition, only: dry_air_composition_update
      use string_utils,    only: int2str
      !------------------------------Arguments----------------------------------
      !(mmr = dry mixing ratio, if not use to_dry_factor to convert)
      real(r8),           intent(in) :: mmr(:,:,:) ! constituents array
      real(r8),           intent(in) :: T(:,:)     ! temperature
      integer,            intent(in) :: lchnk      ! Chunk number
      integer,            intent(in) :: ncol       ! number of columns
      real(r8), optional, intent(in) :: to_dry_factor(:,:)!if mmr moist convert
      !
      !---------------------------Local storage-------------------------------
      real(r8):: sponge_factor(SIZE(mmr, 2))
      character(len=*), parameter :: subname = 'cam_thermo_update: '

      if (present(to_dry_factor)) then
        if (SIZE(to_dry_factor, 1) /= ncol) then
          call endrun(subname//'DIM 1 of to_dry_factor is'//int2str(SIZE(to_dry_factor,1))//'but should be'//int2str(ncol))
        end if
      end if

      sponge_factor = 1.0_r8
      call dry_air_composition_update(mmr, lchnk, ncol, to_dry_factor=to_dry_factor)
      call get_molecular_diff_coef(T(:ncol,:), .true., sponge_factor, kmvis(:ncol,:,lchnk), &
           kmcnd(:ncol,:,lchnk), tracer=mmr(:ncol,:,:), fact=to_dry_factor,  &
           active_species_idx_dycore=thermodynamic_active_species_idx)
    end subroutine cam_thermo_dry_air_update
    !
    !***************************************************************************
    !
    ! cam_thermo_water+update: update water species dependent constants for physics
    !
    !***************************************************************************
    !
    subroutine cam_thermo_water_update(mmr, lchnk, ncol, vcoord, to_dry_factor)
      use air_composition, only: water_composition_update
      !-----------------------------------------------------------------------
      ! Update the physics "constants" that vary
      !-------------------------------------------------------------------------

      !------------------------------Arguments----------------------------------

      real(r8),           intent(in) :: mmr(:,:,:) ! constituents array
      integer,            intent(in) :: lchnk      ! Chunk number
      integer,            intent(in) :: ncol       ! number of columns
      integer,            intent(in) :: vcoord
      real(r8), optional, intent(in) :: to_dry_factor(:,:)
      !
      logical :: lcp

      call water_composition_update(mmr, lchnk, ncol, vcoord, to_dry_factor=to_dry_factor)
    end subroutine cam_thermo_water_update

   !===========================================================================

   !
   !***********************************************************************
   !
   ! Compute enthalpy = cp*T*dp, where dp is pressure level thickness,
   !    cp is generalized cp and T temperature
   !
   ! Note: tracer is in units of m*dp_dry ("mass")
   !
   !***********************************************************************
   !
   subroutine get_enthalpy_1hd(tracer_mass, temp, dp_dry,               &
        enthalpy, active_species_idx_dycore)
      use air_composition, only: dry_air_species_num, get_cp_dry
      ! Dummy arguments
      ! tracer_mass: tracer array (mass weighted)
      real(r8),          intent(in)  :: tracer_mass(:,:,:)
      ! temp: temperature
      real(r8),          intent(in)  :: temp(:,:)
      ! dp_dry: dry presure level thickness
      real(r8),          intent(in)  :: dp_dry(:,:)
      ! enthalpy: enthalpy in each column: sum cp*T*dp
      real(r8),          intent(out) :: enthalpy(:,:)
      !
      ! active_species_idx_dycore:
      !    array of indicies for index of thermodynamic active species in
      !    dycore tracer array (if different from physics index)
      !
      integer, optional, intent(in)  :: active_species_idx_dycore(:)

      ! Local vars
      integer                     :: qdx, itrac
      character(len=*), parameter :: subname = 'get_enthalpy: '

      !
      ! "mass-weighted" cp (dp must be dry)
      !
      if (dry_air_species_num == 0) then
         enthalpy(:,:) = thermodynamic_active_species_cp(0) *         &
              dp_dry(:,:)
      else
         if (present(active_species_idx_dycore)) then
            call get_cp_dry(tracer_mass, active_species_idx_dycore,           &
                 enthalpy, fact=1.0_r8/dp_dry(:,:))
         else
            call get_cp_dry(tracer_mass, thermodynamic_active_species_idx,    &
                 enthalpy, fact=1.0_r8/dp_dry(:,:))
         end if
         enthalpy(:,:) = enthalpy(:,:) * dp_dry(:,:)
      end if
      !
      ! tracer is in units of m*dp ("mass"), where:
      !    m is the dry mixing ratio
      !    dp is the dry pressure level thickness
      !
      do qdx = dry_air_species_num + 1, thermodynamic_active_species_num
         if (present(active_species_idx_dycore)) then
            itrac = active_species_idx_dycore(qdx)
         else
            itrac = thermodynamic_active_species_idx(qdx)
         end if
         enthalpy(:,:) = enthalpy(:,:) +                      &
              (thermodynamic_active_species_cp(qdx) * tracer_mass(:,:,itrac))
      end do
      enthalpy(:,:) = enthalpy(:,:) * temp(:,:)

   end subroutine get_enthalpy_1hd

   !===========================================================================

   subroutine get_enthalpy_2hd(tracer_mass, temp, dp_dry,               &
        enthalpy, active_species_idx_dycore)
      ! Dummy arguments
      ! tracer_mass: tracer array (mass weighted)
      real(r8),          intent(in)  :: tracer_mass(:,:,:,:)
      ! temp: temperature
      real(r8),          intent(in)  :: temp(:,:,:)
      ! dp_dry: dry presure level thickness
      real(r8),          intent(in)  :: dp_dry(:,:,:)
      ! enthalpy: enthalpy in each column: sum cp*T*dp
      real(r8),          intent(out) :: enthalpy(:,:,:)
      !
      ! active_species_idx_dycore:
      !    array of indicies for index of thermodynamic active species in
      !    dycore tracer array (if different from physics index)
      !
      integer, optional, intent(in)  :: active_species_idx_dycore(:)

      ! Local variables
      integer                     :: jdx
      character(len=*), parameter :: subname = 'get_enthalpy_2hd: '

      do jdx = 1, SIZE(tracer_mass, 2)
         call get_enthalpy(tracer_mass(:, jdx, :, :), temp(:, jdx, :),     &
              dp_dry(:, jdx, :), enthalpy(:, jdx, :),                   &
              active_species_idx_dycore=active_species_idx_dycore)
      end do

   end subroutine get_enthalpy_2hd

   !===========================================================================

   !**************************************************************************
   !
   ! get_virtual_temp: Compute virtual temperature T_v
   !
   ! tracer is in units of dry mixing ratio unless optional argument
   !    dp_dry is present in which case tracer is in units of "mass" (=m*dp)
   !
   ! If temperature is not supplied then just return factor that T
   !    needs to be multiplied by to get T_v
   !
   !**************************************************************************
   !
   subroutine get_virtual_temp_1hd(tracer, T_v, temp, dp_dry, sum_q,          &
        active_species_idx_dycore)
      use cam_abortutils,  only: endrun
      use string_utils,    only: int2str
      use air_composition, only: dry_air_species_num, get_R_dry

      ! Dummy Arguments
      ! tracer: tracer array
      real(r8),           intent(in)  :: tracer(:, :, :)
      ! T_v: virtual temperature
      real(r8),           intent(out) :: T_v(:, :)
      ! temp: temperature
      real(r8), optional, intent(in)  :: temp(:, :)
      ! dp_dry: dry pressure level thickness
      real(r8), optional, intent(in)  :: dp_dry(:, :)
      ! sum_q: sum tracer
      real(r8), optional, intent(out) :: sum_q(:, :)
      !
      ! array of indicies for index of thermodynamic active species in
      !    dycore tracer array (if different from physics index)
      !
      integer, optional,  intent(in) :: active_species_idx_dycore(:)

      ! Local Variables
      integer                     :: itrac, qdx
      real(r8)                    :: sum_species(SIZE(tracer, 1), SIZE(tracer, 2))
      real(r8)                    :: factor(SIZE(tracer, 1), SIZE(tracer, 2))
      real(r8)                    :: Rd(SIZE(tracer, 1), SIZE(tracer, 2))
      integer                     :: idx_local(thermodynamic_active_species_num)
      character(len=*), parameter :: subname = 'get_virtual_temp_1hd: '

      if (present(active_species_idx_dycore)) then
         if (SIZE(active_species_idx_dycore) /=                               &
              thermodynamic_active_species_num) then
            call endrun(subname//"SIZE mismatch "//                           &
                 int2str(SIZE(active_species_idx_dycore))//' /= '//           &
                 int2str(thermodynamic_active_species_num))
         end if
         idx_local = active_species_idx_dycore
      else
         idx_local = thermodynamic_active_species_idx
      end if

      call get_sum_species(tracer, idx_local, sum_species, dp_dry=dp_dry, factor=factor)

      call get_R_dry(tracer, idx_local, Rd, fact=factor)
      t_v(:, :) = Rd(:, :)
      do qdx = dry_air_species_num + 1, thermodynamic_active_species_num
         itrac = idx_local(qdx)
         t_v(:, :) = t_v(:, :) + (thermodynamic_active_species_R(qdx) *       &
              tracer(:, :, itrac) * factor(:, :))
      end do
      if (present(temp)) then
         t_v(:, :) = t_v(:, :) * temp(:, :) / (Rd(:, :) * sum_species)
      else
         t_v(:, :) = t_v(:, :) / (Rd(:, :) * sum_species)
      end if
      if (present(sum_q)) then
         sum_q = sum_species
      end if

   end subroutine get_virtual_temp_1hd

   !===========================================================================

   subroutine get_virtual_temp_2hd(tracer, T_v, temp, dp_dry, sum_q,          &
        active_species_idx_dycore)

      ! Dummy Arguments
      ! tracer: tracer array
      real(r8),           intent(in)  :: tracer(:, :, :, :)
      ! T_v: virtual temperature
      real(r8),           intent(out) :: T_v(:, :, :)
      ! temp: temperature
      real(r8), optional, intent(in)  :: temp(:, :, :)
      ! dp_dry: dry pressure level thickness
      real(r8), optional, intent(in)  :: dp_dry(:, :, :)
      ! sum_q: sum tracer
      real(r8), optional, intent(out) :: sum_q(:, :, :)
      !
      ! array of indicies for index of thermodynamic active species in
      !    dycore tracer array (if different from physics index)
      !
      integer, optional,  intent(in) :: active_species_idx_dycore(:)

      ! Local vars
      integer                     :: jdx
      character(len=*), parameter :: subname = 'get_virtual_temp_2hd: '

      ! Rather than do a bunch of copying into temp variables, do the
      !    combinatorics
      do jdx = 1, SIZE(tracer, 2)
         if (present(temp) .and. present(dp_dry) .and. present(sum_q)) then
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 temp=temp(:, jdx, :), dp_dry=dp_dry(:, jdx, :),              &
                 sum_q=sum_q(:, jdx, :),                                      &
                 active_species_idx_dycore=active_species_idx_dycore)
         else if (present(temp) .and. present(dp_dry)) then
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 temp=temp(:, jdx, :), dp_dry=dp_dry(:, jdx, :),              &
                 active_species_idx_dycore=active_species_idx_dycore)
         else if (present(temp) .and. present(sum_q)) then
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 temp=temp(:, jdx, :), sum_q=sum_q(:, jdx, :),                &
                 active_species_idx_dycore=active_species_idx_dycore)
         else if (present(dp_dry) .and. present(sum_q)) then
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 dp_dry=dp_dry(:, jdx, :), sum_q=sum_q(:, jdx, :),            &
                 active_species_idx_dycore=active_species_idx_dycore)
         else if (present(temp)) then
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 temp=temp(:, jdx, :),                                        &
                 active_species_idx_dycore=active_species_idx_dycore)
         else if (present(dp_dry)) then
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 dp_dry=dp_dry(:, jdx, :),                                    &
                 active_species_idx_dycore=active_species_idx_dycore)
         else if (present(sum_q)) then
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 sum_q=sum_q(:, jdx, :),                                      &
                 active_species_idx_dycore=active_species_idx_dycore)
         else
            call get_virtual_temp(tracer(:, jdx, :, :), T_v(:, jdx, :),       &
                 active_species_idx_dycore=active_species_idx_dycore)
         end if
      end do

   end subroutine get_virtual_temp_2hd

   !===========================================================================

   !
   !***************************************************************************
   !
   ! get_sum_species:
   !
   ! Compute sum of thermodynamically active species
   !
   ! tracer is in units of dry mixing ratio unless optional argument
   !    dp_dry is present in which case tracer is in units of "mass" (=m*dp)
   !
   !***************************************************************************
   !
   subroutine get_sum_species_1hd(tracer, active_species_idx,                 &
        sum_species, dp_dry, factor)
      use air_composition, only: dry_air_species_num

      ! Dummy arguments
      ! tracer: Tracer array
      real(r8),           intent(in)  :: tracer(:, :, :)
      ! active_species_idx: Index for thermodynamic active tracers
      integer,            intent(in)  :: active_species_idx(:)
      ! dp_dry: Dry pressure level thickness.
      !         If present, then tracer is in units of mass
      real(r8), optional, intent(in)  :: dp_dry(:, :)
      ! sum_species: sum species
      real(r8),           intent(out) :: sum_species(:, :)
      ! factor: to moist factor 
      real(r8), optional, intent(out) :: factor(:, :)
      ! Local variables
      real(r8) :: factor_loc(SIZE(tracer, 1), SIZE(tracer, 2))
      integer  :: qdx, itrac
      if (present(dp_dry)) then
         factor_loc = 1.0_r8 / dp_dry(:,:)
      else
         factor_loc = 1.0_r8
      end if
      sum_species = 1.0_r8 ! all dry air species sum to 1
      do qdx = dry_air_species_num + 1, thermodynamic_active_species_num
         itrac = active_species_idx(qdx)
         sum_species(:,:) = sum_species(:,:) + (tracer(:,:,itrac) * factor_loc(:,:))
      end do
      if (present(factor)) then
         factor = factor_loc
      end if
   end subroutine get_sum_species_1hd

   !===========================================================================

   subroutine get_sum_species_2hd(tracer, active_species_idx,                 &
        sum_species,dp_dry, factor)

      ! Dummy arguments
      ! tracer: Tracer array
      real(r8),           intent(in)  :: tracer(:, :, :, :)
      ! active_species_idx: Index for thermodynamic active tracers
      integer,            intent(in)  :: active_species_idx(:)
      ! dp_dry: Dry pressure level thickness.
      !         If present, then tracer is in units of mass
      real(r8), optional, intent(in)  :: dp_dry(:, :, :)
      ! sum_species: sum species
      real(r8),           intent(out) :: sum_species(:, :, :)
      ! factor: to moist factor
      real(r8), optional, intent(out) :: factor(:, :, :)
      ! Local variable
      integer                     :: jdx

      do jdx = 1, SIZE(tracer, 2)
         if (present(dp_dry) .and. present(factor)) then
            call get_sum_species(tracer(:, jdx, :, :), active_species_idx,    &
                 sum_species(:, jdx, :), dp_dry=dp_dry(:, jdx, :), factor=factor(:, jdx, :))
         else if (present(dp_dry)) then
            call get_sum_species(tracer(:, jdx, :, :), active_species_idx,    &
                 sum_species(:, jdx, :), dp_dry=dp_dry(:, jdx, :))
         else if (present(factor)) then
            call get_sum_species(tracer(:, jdx, :, :), active_species_idx,    &
                 sum_species(:, jdx, :), factor=factor(:, jdx, :))
         else
            call get_sum_species(tracer(:, jdx, :, :), active_species_idx,    &
                 sum_species(:, jdx, :))
         end if
      end do

   end subroutine get_sum_species_2hd

   !===========================================================================

   !***************************************************************************
   !
   ! get_dp: Compute pressure level thickness from dry pressure and
   !         thermodynamic active species mixing ratios
   !
   ! Tracer can either be in units of dry mixing ratio (mixing_ratio=1) or
   !    "mass" (=m*dp_dry) (mixing_ratio=2)
   !
   !***************************************************************************
   !
   subroutine get_dp_1hd(tracer, mixing_ratio, active_species_idx, dp_dry, dp, ps, ptop)
     use air_composition,  only: dry_air_species_num
     use string_utils,    only: int2str

     real(r8), intent(in)  :: tracer(:, :, :)                    ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is dry mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(:, :)                       ! dry pressure level thickness
     real(r8), intent(out) :: dp(:, :)                           ! pressure level thickness
     real(r8), optional,intent(out) :: ps(:)                     ! surface pressure (if ps present then ptop
                                                                 !                   must be present)
     real(r8), optional,intent(in)  :: ptop                      ! pressure at model top

     integer :: idx, kdx, m_cnst, qdx

     character(len=*), parameter :: subname = 'get_dp_1hd: '

     dp = dp_dry
     if (mixing_ratio == DRY_MIXING_RATIO) then
       do qdx = dry_air_species_num + 1, thermodynamic_active_species_num
         m_cnst = active_species_idx(qdx)
         do kdx = 1, SIZE(tracer, 2)
           do idx = 1, SIZE(tracer, 1)
             dp(idx, kdx) = dp(idx, kdx) + dp_dry(idx, kdx)*tracer(idx, kdx, m_cnst)
           end do
         end do
       end do
     else if (mixing_ratio == MASS_MIXING_RATIO) then
       do qdx = dry_air_species_num + 1, thermodynamic_active_species_num
         m_cnst = active_species_idx(qdx)
         do kdx = 1, SIZE(tracer, 2)
           do idx = 1, SIZE(tracer, 1)
             dp(idx, kdx) = dp(idx, kdx) + tracer(idx, kdx, m_cnst)
           end do
         end do
       end do
     else
        call endrun(subname//'unrecognized input ('//int2str(mixing_ratio)//') for mixing_ratio')
     end if
     if (present(ps)) then
       if (present(ptop)) then
         ps = ptop
         do kdx = 1, SIZE(tracer, 2)
           do idx = 1, SIZE(tracer, 1)
             ps(idx) = ps(idx) + dp(idx, kdx)
           end do
         end do
       else
         call endrun(subname//'if ps is present ptop must be present')
       end if
     end if
   end subroutine get_dp_1hd

   subroutine get_dp_2hd(tracer, mixing_ratio, active_species_idx, dp_dry, dp, ps, ptop)
     ! Version of get_dp for arrays that have a second horizontal index
     real(r8), intent(in)  :: tracer(:,:,:,:)                    ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is dry mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(:,:,:)                      ! dry pressure level thickness
     real(r8), intent(out) :: dp(:,:,:)                          ! pressure level thickness
     real(r8), optional,intent(out) :: ps(:,:)                   ! surface pressure
     real(r8), optional,intent(in)  :: ptop                      ! pressure at model top

     integer :: jdx

     do jdx = 1, SIZE(tracer, 2)
       if (present(ps)) then
          call get_dp(tracer(:, jdx, :, :), mixing_ratio, active_species_idx, &
                  dp_dry(:, jdx, :), dp(:, jdx, :), ps=ps(:,jdx), ptop=ptop)
       else
          call get_dp(tracer(:, jdx, :, :), mixing_ratio, active_species_idx,  &
                  dp_dry(:, jdx, :), dp(:, jdx, :), ptop=ptop)
       end if
     end do

   end subroutine get_dp_2hd
   !===========================================================================

   !*************************************************************************************************************************
   !
   ! compute mid-level (full level) pressure from dry pressure and water tracers
   !
   !*************************************************************************************************************************
   !
   subroutine get_pmid_from_dpdry_1hd(tracer, mixing_ratio, active_species_idx, dp_dry, ptop, pmid, pint, dp)

     real(r8), intent(in)  :: tracer(:,:,:)                      ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(:,:)                        ! dry pressure level thickness
     real(r8), intent(in)  :: ptop                               ! model top pressure
     real(r8), intent(out) :: pmid(:,:)                          ! mid-level pressure
     real(r8), optional, intent(out) :: pint(:,:)                ! half-level pressure
     real(r8), optional, intent(out) :: dp(:,:)                  ! presure level thickness

     real(r8) :: dp_local(SIZE(tracer, 1), SIZE(tracer, 2))      ! local pressure level thickness
     real(r8) :: pint_local(SIZE(tracer, 1), SIZE(tracer, 2) + 1)! local interface pressure

     call get_dp(tracer, mixing_ratio, active_species_idx, dp_dry, dp_local)

     call get_pmid_from_dp(dp_local, ptop, pmid, pint_local)

     if (present(pint)) pint=pint_local
     if (present(dp)) dp=dp_local
   end subroutine get_pmid_from_dpdry_1hd

   !===========================================================================

   !*************************************************************************************************************************
   !
   ! compute mid-level (full level) pressure
   !
   !*************************************************************************************************************************
   !
   subroutine get_pmid_from_dp_1hd(dp, ptop, pmid, pint)
     use dycore, only: dycore_is
     real(r8), intent(in)            :: dp(:,:)     ! pressure level thickness
     real(r8), intent(in)            :: ptop        ! pressure at model top
     real(r8), intent(out)           :: pmid(:,:)   ! mid (full) level pressure
     real(r8), optional, intent(out) :: pint(:,:)   ! pressure at interfaces (half levels)

     real(r8) :: pint_local(SIZE(dp, 1), SIZE(dp,2) + 1)
     integer  :: kdx

     pint_local(:, 1) = ptop
     do kdx = 2, SIZE(dp, 2) + 1
       pint_local(:, kdx) = dp(:, kdx - 1) + pint_local(:, kdx - 1)
     end do

     if (dycore_is('LR') .or. dycore_is('FV3')) then
       do kdx = 1, SIZE(dp, 2)
         pmid(:, kdx) = dp(:, kdx) / (log(pint_local(:, kdx + 1)) - log(pint_local(:, kdx)))
       end do
     else
       do kdx = 1, SIZE(dp, 2)
         pmid(:, kdx) = 0.5_r8 * (pint_local(:, kdx) + pint_local(:, kdx + 1))
       end do
     end if
     if (present(pint)) pint=pint_local
   end subroutine get_pmid_from_dp_1hd

   !===========================================================================

   !****************************************************************************************************************
   !
   ! Compute Exner pressure
   !
   !****************************************************************************************************************
   !
   subroutine get_exner_1hd(tracer, mixing_ratio, active_species_idx, dp_dry, ptop, p00, inv_exner, exner, poverp0)
     use string_utils,    only: int2str
     real(r8), intent(in)  :: tracer(:,:,:)                      ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(:,:)                        ! dry pressure level thickness
     real(r8), intent(in)  :: ptop                               ! pressure at model top
     real(r8), intent(in)  :: p00                                ! reference pressure for Exner pressure (usually 1000hPa)
     logical , intent(in)  :: inv_exner                          ! logical for outputting inverse Exner or Exner pressure
     real(r8), intent(out) :: exner(:,:)
     real(r8), optional, intent(out) :: poverp0(:,:)             ! for efficiency when a routine needs this variable

     real(r8) :: pmid(SIZE(tracer, 1), SIZE(tracer, 2))
     real(r8) :: kappa_dry(SIZE(tracer, 1), SIZE(tracer, 2))
     character(len=*), parameter :: subname = 'get_exner_1hd: '
     !
     ! compute mid level pressure
     !
     call get_pmid_from_dp(tracer, mixing_ratio, active_species_idx, dp_dry, ptop, pmid)
     !
     ! compute kappa = Rd / cpd
     !
     if (mixing_ratio == DRY_MIXING_RATIO) then
       call get_kappa_dry(tracer, active_species_idx, kappa_dry)
     else if (mixing_ratio == MASS_MIXING_RATIO) then
       call get_kappa_dry(tracer, active_species_idx, kappa_dry, 1.0_r8 / dp_dry)
     else
       call endrun(subname//'unrecognized input ('//int2str(mixing_ratio)//') for mixing_ratio')
     end if
     if (inv_exner) then
       exner(:,:) = (p00 / pmid(:,:)) ** kappa_dry(:,:)
     else
       exner(:,:) = (pmid(:,:) / p00) ** kappa_dry(:,:)
     end if
     if (present(poverp0)) poverp0 = pmid(:,:) / p00
   end subroutine get_exner_1hd

   !===========================================================================

   !****************************************************************************************************************
   !
   ! Compute virtual potential temperature from dp_dry, m, T and ptop.
   !
   !****************************************************************************************************************
   !
   subroutine get_virtual_theta_1hd(tracer, mixing_ratio, active_species_idx, dp_dry, ptop, p00, temp, theta_v)
     real(r8), intent(in)  :: tracer(:,:,:)                      ! tracers; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is dry mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(:,:)                        ! dry pressure level thickness
     real(r8), intent(in)  :: ptop                               ! pressure at model top
     real(r8), intent(in)  :: p00                                ! reference pressure for Exner pressure (usually 1000hPa)
     real(r8), intent(in)  :: temp(:,:)                          ! temperature
     real(r8), intent(out) :: theta_v(:,:)                       ! virtual potential temperature

     real(r8) :: iexner(SIZE(tracer, 1), SIZE(tracer, 2))

     call get_exner(tracer, mixing_ratio, active_species_idx, dp_dry, ptop, p00, .true., iexner)

     theta_v(:,:) = temp(:,:) * iexner(:,:)

   end subroutine get_virtual_theta_1hd

   !===========================================================================

   !****************************************************************************************************************
   !
   ! Compute geopotential from dry pressure level thichkness, water tracers, model top pressure and temperature
   !
   !****************************************************************************************************************
   !
   subroutine get_gz_from_dp_dry_ptop_temp_1hd(tracer, mixing_ratio, active_species_idx, &
                     dp_dry, ptop, temp, phis, gz, pmid, dp, T_v)
     use air_composition, only: get_R_dry
     use string_utils,    only: int2str
     real(r8), intent(in)  :: tracer(:,:,:)                      ! tracer; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                       ! 1 => tracer is dry mixing ratio
                                                                 ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)              ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(:,:)                        ! dry pressure level thickness
     real(r8), intent(in)  :: ptop                               ! pressure at model top
     real(r8), intent(in)  :: temp(:,:)                          ! temperature
     real(r8), intent(in)  :: phis(:)                            ! surface geopotential
     real(r8), intent(out) :: gz(:,:)                            ! geopotential
     real(r8), optional, intent(out) :: pmid(:,:)                ! mid-level pressure
     real(r8), optional, intent(out) :: dp(:,:)                  ! pressure level thickness
     real(r8), optional, intent(out) :: t_v(:,:)                 ! virtual temperature


     real(r8), dimension(SIZE(tracer, 1), SIZE(tracer, 2))     :: pmid_local, t_v_local, dp_local, R_dry
     real(r8), dimension(SIZE(tracer, 1), SIZE(tracer, 2) + 1) :: pint
     character(len=*), parameter                               :: subname = 'get_gz_from_dp_dry_ptop_temp_1hd: '
     

     call get_pmid_from_dp(tracer, mixing_ratio, active_species_idx, &
                              dp_dry, ptop, pmid_local, pint=pint, dp=dp_local)
     if (mixing_ratio == DRY_MIXING_RATIO) then
       call get_virtual_temp(tracer, t_v_local, temp=temp, active_species_idx_dycore=active_species_idx)
       call get_R_dry(tracer, active_species_idx, R_dry)
     else if (mixing_ratio == MASS_MIXING_RATIO) then
       call get_virtual_temp(tracer, t_v_local, temp=temp, dp_dry=dp_dry, active_species_idx_dycore=active_species_idx)
       call get_R_dry(tracer,active_species_idx, R_dry, fact=1.0_r8 / dp_dry)
     else
       call endrun(subname//'unrecognized input ('//int2str(mixing_ratio)//') for mixing_ratio')
     end if
     call get_gz(dp_local, T_v_local, R_dry, phis, ptop, gz, pmid_local)

     if (present(pmid)) pmid=pmid_local
     if (present(T_v))  T_v=T_v_local
     if (present(dp))   dp=dp_local
  end subroutine get_gz_from_dp_dry_ptop_temp_1hd

   !===========================================================================

   !***************************************************************************
   !
   ! Compute geopotential from pressure level thickness and virtual temperature
   !
   !***************************************************************************
   !
   subroutine get_gz_given_dp_Tv_Rdry_1hd(dp, T_v, R_dry, phis, ptop, gz, pmid)
     use dycore, only: dycore_is
     real(r8), intent(in)  :: dp   (:,:)                       ! pressure level thickness
     real(r8), intent(in)  :: T_v  (:,:)                       ! virtual temperature
     real(r8), intent(in)  :: R_dry(:,:)                       ! R dry
     real(r8), intent(in)  :: phis (:)                         ! surface geopotential
     real(r8), intent(in)  :: ptop                             ! model top presure
     real(r8), intent(out) :: gz(:,:)                          ! geopotential
     real(r8), optional, intent(out) :: pmid(:,:)              ! mid-level pressure


     real(r8), dimension(SIZE(dp, 1), SIZE(dp, 2))      :: pmid_local
     real(r8), dimension(SIZE(dp, 1), SIZE(dp, 2) + 1)  :: pint
     real(r8), dimension(SIZE(dp, 1))                   :: gzh, Rdry_tv
     integer :: kdx

     call get_pmid_from_dp(dp, ptop, pmid_local, pint)

     !
     ! integrate hydrostatic eqn
     !
     gzh = phis
     if (dycore_is('LR') .or. dycore_is('FV3')) then
       do kdx = SIZE(dp, 2), 1, -1
         Rdry_tv(:) = R_dry(:, kdx) * T_v(:, kdx)
         gz(:, kdx) = gzh(:) + Rdry_tv(:) * (1.0_r8 - pint(:, kdx) / pmid_local(:, kdx))
         gzh(:)  = gzh(:) + Rdry_tv(:) * (log(pint(:, kdx + 1)) - log(pint(:, kdx)))
       end do
     else
       do kdx = SIZE(dp,2), 1, -1
         Rdry_tv(:) = R_dry(:,kdx) * T_v(:, kdx)
         gz(:,kdx) = gzh(:) + Rdry_tv(:) * 0.5_r8 * dp(:, kdx) / pmid_local(:, kdx)
         gzh(:)  = gzh(:) + Rdry_tv(:) * dp(:, kdx) / pmid_local(:, kdx)
       end do
     end if
     if (present(pmid)) pmid=pmid_local
   end subroutine get_gz_given_dp_Tv_Rdry_1hd

   subroutine get_gz_given_dp_Tv_Rdry_2hd(dp, T_v, R_dry, phis, ptop, gz, pmid)
     ! Version of get_gz_given_dp_Tv_Rdry for arrays that have a second horizontal index
     real(r8), intent(in)  :: dp   (:,:,:)                     ! pressure level thickness
     real(r8), intent(in)  :: T_v  (:,:,:)                     ! virtual temperature
     real(r8), intent(in)  :: R_dry(:,:,:)                     ! R dry
     real(r8), intent(in)  :: phis (:,:)                       ! surface geopotential
     real(r8), intent(in)  :: ptop                             ! model top presure
     real(r8), intent(out) :: gz(:,:,:)                        ! geopotential
     real(r8), optional, intent(out) :: pmid(:,:,:)            ! mid-level pressure

     integer :: jdx

     do jdx = 1, SIZE(dp, 2)
       if (present(pmid)) then
         call get_gz(dp(:, jdx, :), T_v(:, jdx, :), R_dry(:, jdx, :), phis(:, jdx), &
              ptop, gz(:, jdx, :), pmid=pmid(:, jdx, :))
       else
         call get_gz(dp(:, jdx, :), T_v(:, jdx, :), R_dry(:, jdx, :), phis(:, jdx), ptop, gz(:, jdx, :))
       end if
     end do


   end subroutine get_gz_given_dp_Tv_Rdry_2hd

   !===========================================================================

   !***************************************************************************
   !
   ! Compute Richardson number at cell interfaces (half levels)
   !
   !***************************************************************************
   !
   subroutine get_Richardson_number_1hd(tracer,mixing_ratio, active_species_idx, dp_dry, ptop, &
        p00, temp, v, Richardson_number, pmid, dp)
     real(r8), intent(in)  :: tracer(:,:,:)                        ! tracer; quantity specified by mixing_ratio arg
     integer,  intent(in)  :: mixing_ratio                         ! 1 => tracer is dry mixing ratio
                                                                   ! 2 => tracer is mass (q*dp)
     integer,  intent(in)  :: active_species_idx(:)                ! index for thermodynamic species in tracer array
     real(r8), intent(in)  :: dp_dry(:,:)                          ! dry pressure level thickness
     real(r8), intent(in)  :: ptop                                 ! pressure at model top
     real(r8), intent(in)  :: p00                                  ! reference pressure for Exner pressure (usually 1000hPa)
     real(r8), intent(in)  :: temp(:,:)                            ! temperature
     real(r8), intent(in)  :: v(:,:,:)                             ! velocity components
     real(r8), intent(out) :: Richardson_number(:,:)
     real(r8), optional, intent(out) :: pmid(:,:)
     real(r8), optional, intent(out) :: dp(:,:)

     real(r8), dimension(SIZE(tracer, 1), SIZE(tracer, 2)) :: gz, theta_v
     real(r8), dimension(SIZE(tracer, 1))                  :: pt1, pt2, phis
     integer :: kdx, kdxm1
     real(r8), parameter:: ustar2 = 1.E-4_r8

     phis = 0.0_r8
     call get_gz(tracer, mixing_ratio, active_species_idx, dp_dry, ptop, temp, phis, gz, pmid=pmid, dp=dp)
     call get_virtual_theta(tracer, mixing_ratio, active_species_idx, dp_dry, ptop, p00, temp, theta_v)
     Richardson_number(:, 1)                   = 0.0_r8
     Richardson_number(:, SIZE(tracer, 2) + 1) = 0.0_r8
     do kdx = SIZE(tracer, 2), 2, -1
       kdxm1 = kdx - 1
       pt1(:) = theta_v(:, kdxm1)
       pt2(:) = theta_v(:, kdx)
       Richardson_number(:, kdx) = (gz(:, kdxm1) - gz(:, kdx)) * (pt1 - pt2) / ( 0.5_r8*(pt1 + pt2) *        &
            ((v(:, 1, kdxm1) - v(:, 1, kdx)) ** 2 + (v(:, 2, kdxm1) - v(:, 2, kdx)) ** 2 + ustar2) )
     end do
   end subroutine get_Richardson_number_1hd

   !
   !****************************************************************************************************************
   !
   ! get surface pressure from dry pressure and thermodynamic active species (e.g., forms of water: water vapor, cldliq, etc.)
   !
   !****************************************************************************************************************
   !
   subroutine get_ps_1hd(tracer_mass, active_species_idx, dp_dry, ps, ptop)
     use air_composition,  only: dry_air_species_num
     
     real(r8), intent(in)   :: tracer_mass(:,:,:)                      ! Tracer array (q*dp)
     real(r8), intent(in)   :: dp_dry(:,:)                             ! dry pressure level thickness
     real(r8), intent(out)  :: ps(:)                                   ! surface pressure
     real(r8), intent(in)   :: ptop
     integer,  intent(in)   :: active_species_idx(:)

     integer                    :: idx, kdx, m_cnst, qdx
     real(r8)                   :: dp(SIZE(tracer_mass, 1), SIZE(tracer_mass, 2))  ! dry pressure level thickness

     dp = dp_dry
     do qdx = dry_air_species_num + 1, thermodynamic_active_species_num
       m_cnst = active_species_idx(qdx)
       do kdx = 1, SIZE(tracer_mass, 2)
         do idx = 1, SIZE(tracer_mass, 1)
           dp(idx, kdx) = dp(idx, kdx) + tracer_mass(idx, kdx, m_cnst)
         end do
       end do
     end do
     ps = ptop
     do kdx = 1, SIZE(tracer_mass, 2)
       do idx = 1, SIZE(tracer_mass, 1)
         ps(idx) = ps(idx) + dp(idx, kdx)
       end do
     end do
   end subroutine get_ps_1hd

   subroutine get_ps_2hd(tracer_mass, active_species_idx, dp_dry, ps, ptop)
     ! Version of get_ps for arrays that have a second horizontal index
     real(r8), intent(in)   :: tracer_mass(:,:,:,:)                      ! Tracer array (q*dp)
     real(r8), intent(in)   :: dp_dry(:,:,:)                             ! dry pressure level thickness
     real(r8), intent(out)  :: ps(:,:)                                   ! surface pressure
     real(r8), intent(in)   :: ptop
     integer,  intent(in)   :: active_species_idx(:)

     integer :: jdx

     do jdx = 1, SIZE(tracer_mass, 2)
       call get_ps(tracer_mass(:, jdx, :, :), active_species_idx, dp_dry(:, jdx, :), ps(:, jdx), ptop)
     end do

   end subroutine get_ps_2hd

   !===========================================================================

   !*************************************************************************************************************************
   !
   ! compute generalized kappa =Rdry/cpdry
   !
   !*************************************************************************************************************************
   !
   subroutine get_kappa_dry_1hd(tracer, active_species_idx, kappa_dry, fact)
     use air_composition,  only: dry_air_species_num, get_R_dry, get_cp_dry
     use physconst,        only: rair, cpair

     real(r8), intent(in)  :: tracer(:,:,:)              !tracer array
     integer,  intent(in)  :: active_species_idx(:)      !index of thermodynamic active tracers
     real(r8), intent(out) :: kappa_dry(:,:)             !kappa dry
     real(r8), optional, intent(in) :: fact(:,:)         !factor for converting tracer to dry mixing ratio
     !
     real(r8), allocatable, dimension(:,:) :: cp_dry,R_dry
     integer                     :: ierr
     character(len=*), parameter :: subname = "get_kappa_dry_1hd"
     character(len=*), parameter :: errstr = subname//": failed to allocate "
     !
     ! dry air not species dependent
     if (dry_air_species_num==0) then
       kappa_dry = rair / cpair
     else
       allocate(R_dry(SIZE(kappa_dry, 1), SIZE(kappa_dry, 2)), stat=ierr)
       if (ierr /= 0) then
         call endrun(errstr//"R_dry")
       end if
       allocate(cp_dry(SIZE(kappa_dry, 1), SIZE(kappa_dry, 2)), stat=ierr)
       if (ierr /= 0) then
         call endrun(errstr//"cp_dry")
       end if
       call get_cp_dry(tracer, active_species_idx, cp_dry, fact=fact)
       call get_R_dry( tracer, active_species_idx, R_dry,  fact=fact)
       kappa_dry = R_dry / cp_dry
       deallocate(R_dry, cp_dry)
     end if
   end subroutine get_kappa_dry_1hd

   subroutine get_kappa_dry_2hd(tracer, active_species_idx, kappa_dry, fact)
     ! Version of get_kappa_dry for arrays that have a second horizontal index
     real(r8), intent(in)  :: tracer(:,:,:,:)              !tracer array
     integer,  intent(in)  :: active_species_idx(:)        !index of thermodynamic active tracers
     real(r8), intent(out) :: kappa_dry(:,:,:)             !kappa dry
     real(r8), optional, intent(in) :: fact(:,:,:)         !factor for converting tracer to dry mixing ratio

     integer :: jdx

     do jdx = 1, SIZE(tracer, 2)
       if (present(fact)) then
         call get_kappa_dry(tracer(:, jdx, :, :), active_species_idx, kappa_dry(:, jdx, :), fact=fact(:, jdx, :))
       else
         call get_kappa_dry(tracer(:, jdx, :, :), active_species_idx, kappa_dry(:, jdx, :))
       end if
     end do

   end subroutine get_kappa_dry_2hd

   !===========================================================================

   !*************************************************************************************************************************
   !
   ! compute reference pressure levels
   !
   !*************************************************************************************************************************
   !
   subroutine get_dp_ref_1hd(hyai, hybi, ps0, phis, dp_ref, ps_ref)
     use physconst,  only: tref, rair
     real(r8), intent(in)           :: hyai(:)
     real(r8), intent(in)           :: hybi(:)
     real(r8), intent(in)           :: ps0
     real(r8), intent(in)           :: phis(:)
     real(r8), intent(out)          :: dp_ref(:,:)
     real(r8), intent(out)          :: ps_ref(:)
     integer :: kdx
     !
     ! use static reference pressure (hydrostatic balance incl. effect of topography)
     !
     ps_ref(:) = ps0 * exp(-phis(:) / (rair * tref))
     do kdx = 1, SIZE(dp_ref, 2)
       dp_ref(:,kdx) = ((hyai(kdx + 1) - hyai(kdx)) * ps0 + (hybi(kdx + 1) - hybi(kdx)) * ps_ref(:))
     end do
   end subroutine get_dp_ref_1hd

   subroutine get_dp_ref_2hd(hyai, hybi, ps0, phis, dp_ref, ps_ref)
     ! Version of get_dp_ref for arrays that have a second horizontal index
     real(r8), intent(in)           :: hyai(:)
     real(r8), intent(in)           :: hybi(:)
     real(r8), intent(in)           :: ps0
     real(r8), intent(in)           :: phis(:,:)
     real(r8), intent(out)          :: dp_ref(:,:,:)
     real(r8), intent(out)          :: ps_ref(:,:)
     integer :: jdx

     do jdx = 1, SIZE(dp_ref, 2)
       call get_dp_ref(hyai, hybi, ps0, phis(:, jdx), dp_ref(:, jdx, :), ps_ref(:, jdx))
     end do

   end subroutine get_dp_ref_2hd

   !===========================================================================

   !*************************************************************************************************************************
   !
   ! compute dry densisty from temperature (temp) and pressure (dp_dry and tracer)
   !
   !*************************************************************************************************************************
   !
   subroutine get_rho_dry_1hd(tracer, temp, ptop, dp_dry, tracer_mass, rho_dry, rhoi_dry, &
              active_species_idx_dycore)
     use air_composition, only: get_R_dry
     ! args
     real(r8), intent(in)           :: tracer(:,:,:)      ! Tracer array
     real(r8), intent(in)           :: temp(:,:)          ! Temperature
     real(r8), intent(in)           :: ptop
     real(r8), intent(in)           :: dp_dry(:,:)
     logical,  intent(in)           :: tracer_mass
     real(r8), optional,intent(out) :: rho_dry(:,:)
     real(r8), optional,intent(out) :: rhoi_dry(:,:)
     !
     ! array of indicies for index of thermodynamic active species in dycore tracer array
     ! (if different from physics index)
     !
     integer, optional, intent(in)   :: active_species_idx_dycore(:)

     ! local vars
     integer :: idx, kdx
     real(r8),  dimension(SIZE(tracer, 1), SIZE(tracer, 2))        :: pmid
     real(r8),  dimension(SIZE(tracer, 1), SIZE(tracer, 2) + 1)    :: pint
     real(r8),  allocatable                                        :: R_dry(:,:)
     integer,  dimension(thermodynamic_active_species_num)         :: idx_local
     integer                     :: ierr
     character(len=*), parameter :: subname = "get_rho_dry_1hd"
     character(len=*), parameter :: errstr = subname//": failed to allocate "

     if (present(active_species_idx_dycore)) then
       idx_local = active_species_idx_dycore
     else
       idx_local = thermodynamic_active_species_idx
     end if
     !
     ! we assume that air is dry where molecular viscosity may be significant
     !
     call get_pmid_from_dp(dp_dry, ptop, pmid, pint=pint)
     if (present(rhoi_dry)) then
       allocate(R_dry(SIZE(tracer, 1), SIZE(tracer, 2) + 1), stat=ierr)
       if (ierr /= 0) then
         call endrun(errstr//"R_dry")
       end if
       if (tracer_mass) then
         call get_R_dry(tracer, idx_local, R_dry, fact=1.0_r8 / dp_dry)
       else
         call get_R_dry(tracer, idx_local, R_dry)
       end if
       do kdx = 2, SIZE(tracer, 2) + 1
         rhoi_dry(:, kdx) = 0.5_r8 * (temp(:, kdx) + temp(:, kdx - 1))!could be more accurate!
         rhoi_dry(:, kdx) = pint(:,kdx) / (rhoi_dry(:, kdx) * R_dry(:, kdx)) !ideal gas law for dry air
       end do
       !
       ! extrapolate top level value
       !
       kdx=1
       rhoi_dry(:, kdx) = 1.5_r8 * (temp(:, kdx) - 0.5_r8 * temp(:, kdx + 1))
       rhoi_dry(:, kdx) = pint(:, kdx) / (rhoi_dry(:, kdx) * R_dry(:, kdx)) !ideal gas law for dry air
       deallocate(R_dry)
     end if
     if (present(rho_dry)) then
       allocate(R_dry(SIZE(tracer, 1), size(rho_dry, 2)), stat=ierr)
       if (ierr /= 0) then
         call endrun(errstr//"R_dry")
       end if
       if (tracer_mass) then
         call get_R_dry(tracer, idx_local, R_dry, fact=1.0_r8 / dp_dry)
       else
         call get_R_dry(tracer, idx_local, R_dry)
       end if
       do kdx = 1, SIZE(rho_dry, 2)
         do idx = 1, SIZE(rho_dry, 1)
           rho_dry(idx, kdx) = pmid(idx, kdx) / (temp(idx, kdx) * R_dry(idx, kdx)) !ideal gas law for dry air
         end do
       end do
       deallocate(R_dry)
     end if
   end subroutine get_rho_dry_1hd

   subroutine get_rho_dry_2hd(tracer, temp, ptop, dp_dry, tracer_mass, rho_dry, rhoi_dry, &
              active_species_idx_dycore)
     ! Version of get_rho_dry for arrays that have a second horizontal index
     real(r8), intent(in)           :: tracer(:,:,:,:)      ! Tracer array
     real(r8), intent(in)           :: temp(:,:,:)          ! Temperature
     real(r8), intent(in)           :: ptop
     real(r8), intent(in)           :: dp_dry(:,:,:)
     logical,  intent(in)           :: tracer_mass
     real(r8), optional,intent(out) :: rho_dry(:,:,:)
     real(r8), optional,intent(out) :: rhoi_dry(:,:,:)
     !
     ! array of indicies for index of thermodynamic active species in dycore tracer array
     ! (if different from physics index)
     !
     integer, optional, intent(in)   :: active_species_idx_dycore(:)

     integer :: jdx

     do jdx = 1, SIZE(tracer, 2)
       if (present(rho_dry) .and. present(rhoi_dry)) then
          call get_rho_dry(tracer(:, jdx, :, :), temp(:, jdx, :), ptop, dp_dry(:, jdx, :), &
              tracer_mass, rho_dry=rho_dry(:, jdx, :), rhoi_dry=rhoi_dry(:, jdx, :), &
              active_species_idx_dycore=active_species_idx_dycore)
       else if (present(rho_dry)) then
          call get_rho_dry(tracer(:, jdx, :, :), temp(:, jdx, :), ptop, dp_dry(:, jdx, :), &
              tracer_mass, rho_dry=rho_dry(:, jdx, :), active_species_idx_dycore=active_species_idx_dycore)
       else if (present(rhoi_dry)) then
          call get_rho_dry(tracer(:, jdx, :, :), temp(:, jdx, :), ptop, dp_dry(:, jdx, :), &
              tracer_mass, rhoi_dry=rhoi_dry(:, jdx, :), active_species_idx_dycore=active_species_idx_dycore)
       else
          call get_rho_dry(tracer(:, jdx, :, :), temp(:, jdx, :), ptop, dp_dry(:, jdx, :), tracer_mass, &
              active_species_idx_dycore=active_species_idx_dycore)
       end if
     end do

   end subroutine get_rho_dry_2hd
   !===========================================================================

   !*************************************************************************************************************************
   !
   ! compute 3D molecular diffusion and thermal conductivity
   !
   !*************************************************************************************************************************
   !
   subroutine get_molecular_diff_coef_1hd(temp, get_at_interfaces, sponge_factor, kmvis, kmcnd, &
        tracer, fact, active_species_idx_dycore, mbarv_in)
     use air_composition,  only: dry_air_species_num, get_mbarv
     use air_composition,  only: kv1, kc1, kv2, kc2, kv_temp_exp, kc_temp_exp

     ! args
     real(r8), intent(in)           :: temp(:,:)                     ! temperature
     logical,  intent(in)           :: get_at_interfaces             ! true: compute kmvis and kmcnd at interfaces
                                                                     ! false: compute kmvis and kmcnd at mid-levels
     real(r8), intent(in)           :: sponge_factor(:)              ! multiply kmvis and kmcnd with sponge_factor
                                                                     ! (for sponge layer)
     real(r8), intent(out)          :: kmvis(:,:)
     real(r8), intent(out)          :: kmcnd(:,:)
     real(r8), intent(in)           :: tracer(:,:,:)                 ! tracer array
     integer,  intent(in), optional :: active_species_idx_dycore(:)  ! index of active species in tracer
     real(r8), intent(in), optional :: fact(:,:)                     ! if tracer is in units of mass or moist
                                                                     ! fact converts to dry mixing ratio: tracer/fact
     real(r8), intent(in), optional :: mbarv_in(:,:)                 ! composition dependent atmosphere mean mass
     !
     ! local vars
     !
     integer :: idx, kdx, icnst, ispecies
     real(r8):: mbarvi, mm, residual             ! Mean mass at mid level
     real(r8):: cnst_vis, cnst_cnd, temp_local
     real(r8), dimension(SIZE(tracer,1), SIZE(sponge_factor, 1)) :: factor, mbarv
     integer,  dimension(thermodynamic_active_species_num)       :: idx_local
     character(len=*), parameter :: subname = 'get_molecular_diff_coef_1hd: '

     !--------------------------------------------
     ! Set constants needed for updates
     !--------------------------------------------

     if (dry_air_species_num==0) then

       cnst_vis = (kv1 * mmro2 * o2_mwi + kv2 * mmrn2 * n2_mwi) * mbar
       cnst_cnd = (kc1 * mmro2 * o2_mwi + kc2 * mmrn2 * n2_mwi) * mbar
       if (get_at_interfaces) then
           do kdx = 2, SIZE(sponge_factor, 1)
             do idx = 1, SIZE(tracer, 1)
               temp_local   = 0.5_r8 * (temp(idx, kdx) + temp(idx, kdx - 1))
               kmvis(idx, kdx) = sponge_factor(kdx) * cnst_vis * temp_local ** kv_temp_exp
               kmcnd(idx, kdx) = sponge_factor(kdx) * cnst_cnd * temp_local ** kc_temp_exp
             end do
           end do
           !
           ! extrapolate top level value
           !
           kmvis(1:SIZE(tracer, 1), 1) = 1.5_r8 * kmvis(1:SIZE(tracer, 1), 2) - 0.5_r8 * kmvis(1:SIZE(tracer, 1), 3)
           kmcnd(1:SIZE(tracer, 1), 1) = 1.5_r8 * kmcnd(1:SIZE(tracer, 1), 2) - 0.5_r8 * kmcnd(1:SIZE(tracer, 1), 3)
       else if (.not. get_at_interfaces) then
         do kdx = 1, SIZE(sponge_factor, 1)
           do idx = 1, SIZE(tracer, 1)
             kmvis(idx, kdx) = sponge_factor(kdx) * cnst_vis * temp(idx, kdx) ** kv_temp_exp
             kmcnd(idx, kdx) = sponge_factor(kdx) * cnst_cnd * temp(idx, kdx) ** kc_temp_exp
           end do
         end do
       else
         call endrun(subname//'get_at_interfaces must be .true. or .false.')
       end if
     else
       if (present(active_species_idx_dycore)) then
         idx_local = active_species_idx_dycore
       else
         idx_local = thermodynamic_active_species_idx
       end if
       if (present(fact)) then
         factor = fact(:,:)
       else
         factor = 1.0_r8
       endif
       if (present(mbarv_in)) then
         mbarv = mbarv_in
       else
         call get_mbarv(tracer, idx_local, mbarv, fact=factor)
       end if
       !
       ! major species dependent code
       !
       if (get_at_interfaces) then
         do kdx = 2, SIZE(sponge_factor, 1)
           do idx = 1, SIZE(tracer, 1)
             kmvis(idx, kdx) = 0.0_r8
             kmcnd(idx, kdx) = 0.0_r8
             residual = 1.0_r8
             do icnst = 1, dry_air_species_num
               ispecies = idx_local(icnst)
               mm       = 0.5_r8 * (tracer(idx, kdx, ispecies) * factor(idx, kdx) + &
                          tracer(idx, kdx - 1, ispecies) * factor(idx, kdx-1))
               kmvis(idx, kdx) = kmvis(idx, kdx) + thermodynamic_active_species_kv(icnst) * &
                                           thermodynamic_active_species_mwi(icnst) * mm
               kmcnd(idx, kdx) = kmcnd(idx, kdx) + thermodynamic_active_species_kc(icnst) * &
                                           thermodynamic_active_species_mwi(icnst) * mm
               residual        = residual - mm
             end do
             icnst = 0 ! N2
             kmvis(idx, kdx) = kmvis(idx, kdx) + thermodynamic_active_species_kv(icnst) * &
                                         thermodynamic_active_species_mwi(icnst) * residual
             kmcnd(idx, kdx) = kmcnd(idx, kdx) + thermodynamic_active_species_kc(icnst) * &
                                         thermodynamic_active_species_mwi(icnst) * residual

             temp_local = 0.5_r8 * (temp(idx, kdx - 1) + temp(idx, kdx))
             mbarvi = 0.5_r8 * (mbarv(idx, kdx - 1) + mbarv(idx, kdx))
             kmvis(idx, kdx) = kmvis(idx, kdx) * mbarvi * temp_local ** kv_temp_exp
             kmcnd(idx, kdx) = kmcnd(idx, kdx) * mbarvi * temp_local ** kc_temp_exp
           enddo
         end do
         do idx = 1, SIZE(tracer, 1)
           kmvis(idx, 1) = 1.5_r8 * kmvis(idx, 2) - .5_r8 * kmvis(idx, 3)
           kmcnd(idx, 1) = 1.5_r8 * kmcnd(idx, 2) - .5_r8 * kmcnd(idx, 3)
           kmvis(idx, SIZE(sponge_factor, 1) + 1) = kmvis(idx, SIZE(sponge_factor, 1))
           kmcnd(idx, SIZE(sponge_factor, 1) + 1) = kmcnd(idx, SIZE(sponge_factor, 1))
         end do
       else if (.not. get_at_interfaces) then
         do kdx = 1, SIZE(sponge_factor, 1)
           do idx = 1, SIZE(tracer, 1)
             kmvis(idx, kdx) = 0.0_r8
             kmcnd(idx, kdx) = 0.0_r8
             residual = 1.0_r8
             do icnst = 1, dry_air_species_num - 1
               ispecies = idx_local(icnst)
               mm = tracer(idx, kdx, ispecies) * factor(idx, kdx)
               kmvis(idx, kdx) = kmvis(idx, kdx) + thermodynamic_active_species_kv(icnst) * &
                                 thermodynamic_active_species_mwi(icnst) * mm
               kmcnd(idx, kdx) = kmcnd(idx, kdx) + thermodynamic_active_species_kc(icnst) * &
                                 thermodynamic_active_species_mwi(icnst) * mm
               residual        = residual - mm
             end do
             icnst = dry_air_species_num
             kmvis(idx, kdx) = kmvis(idx, kdx) + thermodynamic_active_species_kv(icnst) * &
                               thermodynamic_active_species_mwi(icnst) * residual
             kmcnd(idx, kdx) = kmcnd(idx, kdx) + thermodynamic_active_species_kc(icnst) * &
                               thermodynamic_active_species_mwi(icnst) * residual

             kmvis(idx, kdx) = kmvis(idx, kdx) * mbarv(idx, kdx) * temp(idx, kdx) ** kv_temp_exp
             kmcnd(idx, kdx) = kmcnd(idx, kdx) * mbarv(idx, kdx) * temp(idx, kdx) ** kc_temp_exp
           end do
         end do
       else
         call endrun(subname//'get_at_interfaces must be .true. or .false.')
       end if
     end if
   end subroutine get_molecular_diff_coef_1hd

   subroutine get_molecular_diff_coef_2hd(temp, get_at_interfaces, sponge_factor, kmvis, kmcnd, &
        tracer, fact, active_species_idx_dycore, mbarv_in)
     ! Version of get_molecular_diff_coef for arrays that have a second horizontal index
     real(r8), intent(in)           :: temp(:,:,:)                     ! temperature
     logical,  intent(in)           :: get_at_interfaces               ! true: compute kmvis and kmcnd at interfaces
                                                                       ! false: compute kmvis and kmcnd at mid-levels
     real(r8), intent(in)           :: sponge_factor(:)                ! multiply kmvis and kmcnd with sponge_factor
                                                                       ! (for sponge layer)
     real(r8), intent(out)          :: kmvis(:,:,:)
     real(r8), intent(out)          :: kmcnd(:,:,:)
     real(r8), intent(in)           :: tracer(:,:,:,:)                 ! tracer array
     integer,  intent(in), optional :: active_species_idx_dycore(:)    ! index of active species in tracer
     real(r8), intent(in), optional :: fact(:,:,:)                     ! if tracer is in units of mass or moist
                                                                       ! fact converts to dry mixing ratio: tracer/fact
     real(r8), intent(in), optional :: mbarv_in(:,:,:)                 ! composition dependent atmosphere mean mass
     integer :: jdx

     do jdx = 1, SIZE(tracer, 2)
       if (present(fact) .and. present(mbarv_in)) then
         call get_molecular_diff_coef(temp(:, jdx, :), get_at_interfaces, sponge_factor, &
              kmvis(:, jdx, :), kmcnd(:, jdx, :), tracer(:, jdx, :, :), fact=fact(:, jdx, :), &
              active_species_idx_dycore=active_species_idx_dycore, mbarv_in=mbarv_in(:, jdx, :))
       else if (present(fact)) then
         call get_molecular_diff_coef(temp(:, jdx, :), get_at_interfaces, sponge_factor, &
              kmvis(:, jdx, :), kmcnd(:, jdx, :), tracer(:, jdx, :, :), fact=fact(:, jdx, :), &
              active_species_idx_dycore=active_species_idx_dycore)
       else if (present(mbarv_in)) then
         call get_molecular_diff_coef(temp(:, jdx, :), get_at_interfaces, sponge_factor, &
              kmvis(:, jdx, :), kmcnd(:, jdx, :), tracer(:, jdx, :, :), &
              active_species_idx_dycore=active_species_idx_dycore, mbarv_in=mbarv_in(:, jdx, :))
       else
         call get_molecular_diff_coef(temp(:, jdx, :), get_at_interfaces, sponge_factor, &
              kmvis(:, jdx, :), kmcnd(:, jdx, :), tracer(:, jdx, :, :), &
              active_species_idx_dycore=active_species_idx_dycore)
       end if
     end do

   end subroutine get_molecular_diff_coef_2hd
   !===========================================================================

   !***************************************************************************
   !
   ! compute reference vertical profile of density, molecular diffusion and thermal conductivity
   !
   !***************************************************************************
   !
   subroutine get_molecular_diff_coef_reference(tref,press,sponge_factor,kmvis_ref,kmcnd_ref,rho_ref)
     use physconst,       only: rair
     use air_composition, only: kv1, kv2, kc1, kc2, kv_temp_exp, kc_temp_exp
     ! args
     real(r8), intent(in)           :: tref                 !reference temperature
     real(r8), intent(in)           :: press(:)             !pressure
     real(r8), intent(in)           :: sponge_factor(:)     !multiply kmvis and kmcnd with sponge_factor (for sponge layer)
     real(r8), intent(out)          :: kmvis_ref(:)         !reference molecular diffusion coefficient
     real(r8), intent(out)          :: kmcnd_ref(:)         !reference thermal conductivity coefficient
     real(r8), intent(out)          :: rho_ref(:)           !reference density

     ! local vars
     integer :: kdx

     !--------------------------------------------
     ! Set constants needed for updates
     !--------------------------------------------

     do kdx = 1, SIZE(press, 1)
       rho_ref(kdx) = press(kdx) / (tref * rair) !ideal gas law for dry air
       kmvis_ref(kdx) = sponge_factor(kdx) * &
            (kv1 * mmro2 * o2_mwi +         &
             kv2 * mmrn2 * n2_mwi) * mbar *    &
             tref ** kv_temp_exp
       kmcnd_ref(kdx) = sponge_factor(kdx) * &
            (kc1 * mmro2 * o2_mwi +         &
             kc2 * mmrn2 * n2_mwi) * mbar *    &
             tref ** kc_temp_exp
     end do
   end subroutine get_molecular_diff_coef_reference

   !==========================================================================

   !
   !***************************************************************************
   !
   ! cam_thermo_calc_kappav: update species dependent kappa for FV dycore
   !
   !***************************************************************************
   !
   subroutine cam_thermo_calc_kappav_2hd(tracer, kappav, cpv)
      use air_composition, only: get_R_dry, get_cp_dry
      ! assumes moist MMRs

      ! Dummy arguments
      real(r8),           intent(in)  :: tracer(:, :, :, :)
      real(r8),           intent(out) :: kappav(:, :, :)
      real(r8), optional, intent(out) :: cpv(:, :, :)

      ! Local variables
      real(r8) :: rgas_var(SIZE(tracer, 1), SIZE(tracer, 2), SIZE(tracer, 3))
      real(r8) :: cp_var(SIZE(tracer, 1), SIZE(tracer, 2), SIZE(tracer, 3))
      integer :: ind, jnd, knd

      !-----------------------------------------------------------------------
      !  Calculate constituent dependent specific heat, gas constant and cappa
      !-----------------------------------------------------------------------
      call get_R_dry(tracer, thermodynamic_active_species_idx, rgas_var)
      call get_cp_dry(tracer, thermodynamic_active_species_idx, cp_var)
      !$omp parallel do private(ind,jnd,knd)
      do knd = 1, SIZE(tracer, 3)
         do jnd = 1, SIZE(tracer, 2)
            do ind = 1, SIZE(tracer, 1)
               kappav(ind,jnd,knd) = rgas_var(ind,jnd,knd) / cp_var(ind,jnd,knd)
            end do
         end do
      end do

      if (present(cpv)) then
         cpv(:,:,:) = cp_var(:,:,:)
      end if

   end subroutine cam_thermo_calc_kappav_2hd

   !===========================================================================
   !
   !***************************************************************************
   !
   ! compute column integrated total energy consistent with vertical
   !    coordinate as well as vertical integrals of water mass (H2O,wv,liq,ice)
   !
   ! if subroutine is asked to compute "te" then the latent heat terms are
   !    added to the kinetic (ke), internal + geopotential (se)  energy terms
   !
   ! subroutine assumes that enthalpy term (rho*cp*T) uses dry air heat capacity
   !
   !***************************************************************************
   !
   subroutine get_hydrostatic_energy_1hd(tracer, moist_mixing_ratio, pdel_in, &
        cp_or_cv, U, V, T, vcoord, ptop, phis, z_mid, dycore_idx, qidx,       &
        te, se, po, ke, wv, H2O, liq, ice)

      use cam_logfile,     only: iulog
      use dyn_tests_utils, only: vc_height, vc_moist_pressure, vc_dry_pressure
      use air_composition, only: wv_idx
      use physconst,       only: rga, latvap, latice

      ! Dummy arguments
      ! tracer: tracer mixing ratio
      !
      ! note - if pdeldry passed to subroutine then tracer mixing ratio must be dry
      real(r8), intent(in)            :: tracer(:,:,:)
      logical, intent(in)             :: moist_mixing_ratio
      ! pdel: pressure level thickness
      real(r8), intent(in)            :: pdel_in(:,:)
      ! cp_or_cv: dry air heat capacity under constant pressure or
      !           constant volume (depends on vcoord)
      real(r8), intent(in)            :: cp_or_cv(:,:)
      real(r8), intent(in)            :: U(:,:)
      real(r8), intent(in)            :: V(:,:)
      real(r8), intent(in)            :: T(:,:)
      integer,  intent(in)            :: vcoord ! vertical coordinate
      real(r8), intent(in),  optional :: ptop(:)
      real(r8), intent(in),  optional :: phis(:)
      real(r8), intent(in),  optional :: z_mid(:,:)
      ! dycore_idx: use dycore index for thermodynamic active species
      logical,  intent(in),  optional :: dycore_idx
      ! qidx: Index of water vapor
      integer,  intent(in),  optional :: qidx
      ! H2O: vertically integrated total water
      real(r8), intent(out), optional :: H2O(:)
      ! TE: vertically integrated total energy
      real(r8), intent(out), optional :: te (:)
      ! KE: vertically integrated kinetic energy
      real(r8), intent(out), optional :: ke (:)
      ! SE: vertically integrated enthalpy (pressure coordinate) 
      !     or internal energy (z coordinate)
      real(r8), intent(out), optional :: se (:)
      ! PO: vertically integrated PHIS term (pressure coordinate)
      !     or potential energy (z coordinate)
      real(r8), intent(out), optional :: po (:)
      ! WV: vertically integrated water vapor
      real(r8), intent(out), optional :: wv (:)
      ! liq: vertically integrated liquid
      real(r8), intent(out), optional :: liq(:)
      ! ice: vertically integrated ice
      real(r8), intent(out), optional :: ice(:)

      ! Local variables
      real(r8) :: ke_vint(SIZE(tracer, 1))  ! Vertical integral of KE
      real(r8) :: se_vint(SIZE(tracer, 1))  ! Vertical integral of enthalpy or internal energy
      real(r8) :: po_vint(SIZE(tracer, 1))  ! Vertical integral of PHIS or potential energy
      real(r8) :: wv_vint(SIZE(tracer, 1))  ! Vertical integral of wv
      real(r8) :: liq_vint(SIZE(tracer, 1)) ! Vertical integral of liq
      real(r8) :: ice_vint(SIZE(tracer, 1)) ! Vertical integral of ice
      real(r8) :: pdel(SIZE(tracer, 1),SIZE(tracer, 2)) !moist pressure level thickness
      real(r8)                      :: latsub ! latent heat of sublimation

      integer                       :: ierr
      integer                       :: kdx, idx ! coord indices
      integer                       :: qdx      ! tracer index
      integer                       :: wvidx    ! water vapor index
      integer,          allocatable :: species_idx(:)
      integer,          allocatable :: species_liq_idx(:)
      integer,          allocatable :: species_ice_idx(:)
      character(len=*), parameter   :: subname = 'get_hydrostatic_energy'

      allocate(species_idx(thermodynamic_active_species_num), stat=ierr)
      if ( ierr /= 0 ) then
         call endrun(subname//': allocation error for species_idx array')
      end if
      allocate(species_liq_idx(thermodynamic_active_species_liq_num), stat=ierr)
      if ( ierr /= 0 ) then
         call endrun(subname//': allocation error for species_liq_idx array')
      end if
      allocate(species_ice_idx(thermodynamic_active_species_ice_num), stat=ierr)
      if ( ierr /= 0 ) then
         call endrun(subname//': allocation error for species_ice_idx array')
      end if

      if (present(dycore_idx))then
         if (dycore_idx) then
            species_idx(:) = thermodynamic_active_species_idx_dycore(:)
            species_liq_idx(:) = thermodynamic_active_species_liq_idx_dycore(:)
            species_ice_idx(:) = thermodynamic_active_species_ice_idx_dycore(:)
         else
            species_idx(:) = thermodynamic_active_species_idx(:)
            species_liq_idx(:) = thermodynamic_active_species_liq_idx(:)
            species_ice_idx(:) = thermodynamic_active_species_ice_idx(:)
         end if
      else
         species_idx(:) = thermodynamic_active_species_idx(:)
         species_liq_idx(:) = thermodynamic_active_species_liq_idx(:)
         species_ice_idx(:) = thermodynamic_active_species_ice_idx(:)
      end if

      if (present(qidx)) then
         wvidx = qidx
      else
         wvidx = wv_idx
      end if

      if (moist_mixing_ratio) then
        pdel     = pdel_in
      else
        pdel     = pdel_in
        do qdx = dry_air_species_num+1, thermodynamic_active_species_num
          pdel(:,:) = pdel(:,:) + pdel_in(:, :)*tracer(:,:,species_idx(qdx))
        end do
      end if

      ke_vint = 0._r8
      se_vint = 0._r8
      select case (vcoord)
      case(vc_moist_pressure, vc_dry_pressure)
         if (.not. present(ptop).or. (.not. present(phis))) then
            write(iulog, *) subname, ' ptop and phis must be present for ',     &
                 'moist/dry pressure vertical coordinate'
            call endrun(subname//':  ptop and phis must be present for '//      &
                 'moist/dry pressure vertical coordinate')
         end if
         po_vint = ptop
         do kdx = 1, SIZE(tracer, 2)
            do idx = 1, SIZE(tracer, 1)
               ke_vint(idx) = ke_vint(idx) + (pdel(idx, kdx) *                &
                    0.5_r8 * (U(idx, kdx)**2 + V(idx, kdx)**2)) * rga
               se_vint(idx) = se_vint(idx) + (T(idx, kdx) *                   &
                    cp_or_cv(idx, kdx) * pdel(idx, kdx) * rga)
               po_vint(idx) =  po_vint(idx)+pdel(idx, kdx)

            end do
         end do
         do idx = 1, SIZE(tracer, 1)
            po_vint(idx) =  (phis(idx) * po_vint(idx) * rga)
         end do
      case(vc_height)
         if (.not. present(phis)) then
            write(iulog, *) subname, ' phis must be present for ',     &
                 'heigt-based vertical coordinate'
            call endrun(subname//':  phis must be present for '//      &
                 'height-based vertical coordinate')
         end if
         po_vint = 0._r8
         do kdx = 1, SIZE(tracer, 2)
            do idx = 1, SIZE(tracer, 1)
               ke_vint(idx) = ke_vint(idx) + (pdel(idx, kdx) *                &
                    0.5_r8 * (U(idx, kdx)**2 + V(idx, kdx)**2) * rga)
               se_vint(idx) = se_vint(idx) + (T(idx, kdx) *                   &
                    cp_or_cv(idx, kdx) * pdel(idx, kdx) * rga)
               ! z_mid is height above ground
               po_vint(idx) = po_vint(idx) + (z_mid(idx, kdx) +               &
                    phis(idx) * rga) * pdel(idx, kdx)
            end do
         end do
      case default
         write(iulog, *) subname, ' vertical coordinate not supported: ', vcoord
         call endrun(subname//': vertical coordinate not supported')
      end select
      if (present(te)) then
         te  = se_vint + po_vint+ ke_vint
      end if
      if (present(se)) then
         se = se_vint
      end if
      if (present(po)) then
         po = po_vint
      end if
      if (present(ke)) then
         ke = ke_vint
      end if
      !
      ! vertical integral of total liquid water
      !
      if (.not.moist_mixing_ratio) then
        pdel = pdel_in! set pseudo density to dry
      end if

      wv_vint = 0._r8
      do kdx = 1, SIZE(tracer, 2)
        do idx = 1, SIZE(tracer, 1)
          wv_vint(idx) = wv_vint(idx) + (tracer(idx, kdx, wvidx) *       &
               pdel(idx, kdx) * rga)
        end do
      end do
      if (present(wv)) wv = wv_vint

      liq_vint = 0._r8
      do qdx = 1, thermodynamic_active_species_liq_num
         do kdx = 1, SIZE(tracer, 2)
            do idx = 1, SIZE(tracer, 1)
               liq_vint(idx) = liq_vint(idx) + (pdel(idx, kdx) *         &
                    tracer(idx, kdx, species_liq_idx(qdx)) * rga)
            end do
         end do
      end do
      if (present(liq)) liq = liq_vint

      !
      ! vertical integral of total frozen (ice) water
      !
      ice_vint = 0._r8
      do qdx = 1, thermodynamic_active_species_ice_num
         do kdx = 1, SIZE(tracer, 2)
            do idx = 1, SIZE(tracer, 1)
               ice_vint(idx) = ice_vint(idx) + (pdel(idx, kdx) *              &
                    tracer(idx, kdx, species_ice_idx(qdx))  * rga)
            end do
         end do
      end do
      if (present(ice)) ice = ice_vint
      ! Compute vertical integrals of total water.
      if (present(H2O)) then
         H2O = wv_vint + liq_vint + ice_vint
      end if
      !
      ! latent heat terms depend on enthalpy reference state
      !
      latsub = latvap + latice
      if (present(te)) then
         select case (TRIM(enthalpy_reference_state))
         case('ice')
            te = te + (latsub * wv_vint) + (latice * liq_vint)
         case('liq')
            te = te + (latvap * wv_vint) - (latice * ice_vint)
         case('wv')
            te = te - (latvap * liq_vint) - (latsub * ice_vint)
         case default
            write(iulog, *) subname, ' enthalpy reference state not ',        &
                 'supported: ', TRIM(enthalpy_reference_state)
            call endrun(subname//': enthalpy reference state not supported')
         end select
      end if
      deallocate(species_idx, species_liq_idx, species_ice_idx)
    end subroutine get_hydrostatic_energy_1hd

end module cam_thermo
