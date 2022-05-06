module cam_thermo

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use air_composition, only: thermodynamic_active_species_num
   use air_composition, only: thermodynamic_active_species_idx
   use air_composition, only: thermodynamic_active_species_idx_dycore
   use air_composition, only: thermodynamic_active_species_cp
   use air_composition, only: thermodynamic_active_species_cv
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
   use air_composition, only: enthalpy_reference_state
   use air_composition, only: mmro2, mmrn2, o2_mwi, n2_mwi, mbar

   implicit none
   private
   save

   ! subroutines to compute thermodynamic quantities
   !
   ! See Lauritzen et al. (2018) for formulaes
   !     DOI: 10.1029/2017MS001257
   !     https://opensky.ucar.edu/islandora/object/articles:21929

   ! cam_thermo_init: Initialize constituent dependent properties
   public :: cam_thermo_init
   ! cam_thermo_update: Update constituent dependent properties
   public :: cam_thermo_update
   ! get_cp_dry: (generalized) heat capacity for dry air
   public :: get_cp_dry
   ! get_cp: (generalized) heat capacity
   public :: get_cp
   ! get_R_dry: (generalized) dry air gas constant
   public :: get_R_dry
   ! get_thermal_energy: thermal energy quantity = dp*cp*T
   public :: get_thermal_energy
   ! get_virtual_temp: virtual temperature
   public :: get_virtual_temp
   ! get_sum_species: sum of thermodynamically active species:
   !                  Note: dp = dp_dry * sum_species
   public :: get_sum_species
   ! get_virtual_theta: virtual potential temperature
   public :: get_virtual_theta
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
   ! get_rho_dry: dry densisty from temperature (temp) and
   !              pressure (dp_dry and tracer)
   public :: get_rho_dry
   ! get_exner: Exner pressure
   public :: get_exner
   ! get_hydrostatic_energy: Vertically integrated total energy
   public :: get_hydrostatic_energy

   ! Public variables
   !---------------  Variables below here are for WACCM-X ---------------------
   ! cpairv:  composition dependent specific heat at constant pressure
   real(r8), public, protected, allocatable :: cpairv(:,:,:)
   ! rairv: composition dependent gas "constant"
   real(r8), public, protected, allocatable :: rairv(:,:,:)
   ! cappav: rairv / cpairv
   real(r8), public, protected, allocatable :: cappav(:,:,:)
   ! mbarv: composition dependent atmosphere mean mass
   real(r8), public, protected, allocatable :: mbarv(:,:,:)
   ! kmvis: molecular viscosity      kg/m/s
   real(r8), public, protected, allocatable :: kmvis(:,:,:)
   ! kmcnd: molecular conductivity   J/m/s/K
   real(r8), public, protected, allocatable :: kmcnd(:,:,:)

   !------------- Variables for consistent themodynamics --------------------
   !

   ! coefficients in expressions for molecular diffusion coefficients
   ! kv1,..,kv4 are coefficients for kmvis calculation
   ! kc1,..,kc4 are coefficients for kmcnd calculation
   real(r8), parameter :: kv1 = 4.03_r8
   real(r8), parameter :: kv2 = 3.42_r8
   real(r8), parameter :: kv3 = 3.9_r8
   real(r8), parameter :: kv4 = 0.69_r8
   real(r8), parameter :: kc1 = 56._r8
   real(r8), parameter :: kc2 = 56._r8
   real(r8), parameter :: kc3 = 75.9_r8
   real(r8), parameter :: kc4 = 0.69_r8

   !
   ! Interfaces for public routines
   interface get_cp_dry
      module procedure get_cp_dry_1hd
      module procedure get_cp_dry_2hd
   end interface get_cp_dry

   interface get_cp
      module procedure get_cp_1hd
      module procedure get_cp_2hd
   end interface get_cp

   interface get_R_dry
      module procedure get_R_dry_1hd
      module procedure get_R_dry_2hd
   end interface get_R_dry

   interface get_R
      module procedure get_R_1hd
      module procedure get_R_2hd
   end interface get_R

   interface get_gz
      ! get_gz_geopotential (with dp_dry, ptop, temp, and phis as input)
      module procedure get_gz_from_dp_dry_ptop_temp
      ! get_gz_given_dp_Tv_Rdry: geopotential (with dp,dry R and Tv as input)
      module procedure get_gz_given_dp_Tv_Rdry
   end interface get_gz

   interface get_thermal_energy
      module procedure get_thermal_energy_1hd
      module procedure get_thermal_energy_2hd
   end interface get_thermal_energy

   interface get_virtual_temp
      module procedure get_virtual_temp_1hd
      module procedure get_virtual_temp_2hd
   end interface get_virtual_temp

   interface get_sum_species
      module procedure get_sum_species_1hd
      module procedure get_sum_species_2hd
   end interface get_sum_species

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

CONTAINS

   !===========================================================================

   subroutine cam_thermo_init()
      use shr_infnan_mod,  only: assignment(=), shr_infnan_qnan
      use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
      use physconst,       only: cpair, rair, mwdry

      integer                     :: ierr
      character(len=*), parameter :: subname = "cam_thermo_init"
      character(len=*), parameter :: errstr = subname//": failed to allocate "

      !------------------------------------------------------------------------
      !  Allocate constituent dependent properties
      !------------------------------------------------------------------------
      allocate(cpairv(pcols,pver,begchunk:endchunk), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"cpairv")
      end if
      allocate(rairv(pcols,pver,begchunk:endchunk),  stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"rairv")
      end if
      allocate(cappav(pcols,pver,begchunk:endchunk), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"cappav")
      end if
      allocate(mbarv(pcols,pver,begchunk:endchunk),  stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"mbarv")
      end if
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
      cpairv(:pcols, :pver, begchunk:endchunk) = cpair
      rairv(:pcols,  :pver, begchunk:endchunk) = rair
      cappav(:pcols, :pver, begchunk:endchunk) = rair / cpair
      mbarv(:pcols,  :pver, begchunk:endchunk) = mwdry
      kmvis(:pcols,  :pver, begchunk:endchunk) = shr_infnan_qnan
      kmcnd(:pcols,  :pver, begchunk:endchunk) = shr_infnan_qnan

   end subroutine cam_thermo_init

   !===========================================================================

   !***************************************************************************
   !
   ! cam_thermo_update: update species dependent constants for physics
   !
   !***************************************************************************
   !
   subroutine cam_thermo_update(mmr, T, lchnk, ncol, to_moist_factor)

      !-----------------------------------------------------------------------
      ! Update the physics "constants" that vary
      !-------------------------------------------------------------------------

      !------------------------------Arguments----------------------------------

      real(r8),           intent(in) :: mmr(:,:,:) ! constituents array
      real(r8),           intent(in) :: T(:,:)     ! temperature
      integer,            intent(in) :: lchnk      ! Chunk number
      integer,            intent(in) :: ncol       ! number of columns
      real(r8), optional, intent(in) :: to_moist_factor(:,:)
      !
      !---------------------------Local storage-------------------------------
      real(r8):: to_moist_fact(SIZE(mmr, 1), SIZE(mmr, 2))
      real(r8):: sponge_factor(SIZE(mmr, 2))


      if (present(to_moist_factor)) then
         to_moist_fact(:ncol,:) = to_moist_factor(:ncol,:)
      else
         to_moist_fact(:,:) = 1._r8
      end if

      !--------------------------------------------
      ! update cpairv, rairv, mbarv, and cappav
      !--------------------------------------------
      call get_R_dry(mmr(:ncol,:,:), thermodynamic_active_species_idx, &
           rairv(:ncol,:,lchnk), fact=to_moist_fact(:ncol,:))
      call get_cp_dry(mmr(:ncol,:,:), thermodynamic_active_species_idx, &
           cpairv(:ncol,:,lchnk), fact=to_moist_fact(:ncol,:))
      call get_mbarv(mmr(:ncol,:,:), thermodynamic_active_species_idx, &
           mbarv(:ncol,:,lchnk), fact=to_moist_fact(:ncol,:))
      sponge_factor = 1.0_r8
      call get_molecular_diff_coef(T(:ncol,:), 1, sponge_factor,kmvis(:ncol,:,lchnk), &
           kmcnd(:ncol,:,lchnk), pcnst, tracer=mmr(:ncol,:,:), fact=to_moist_fact(:ncol,:),                &
           active_species_idx_dycore=thermodynamic_active_species_idx)
      cappav(:ncol,:,lchnk) = rairv(:ncol,:,lchnk) / cpairv(:ncol,:,lchnk)

   end subroutine cam_thermo_update

   !===========================================================================

   !***************************************************************************
   !
   ! get_cp_dry: Compute dry air heat capacity under constant pressure
   !
   !***************************************************************************
   !
   subroutine get_cp_dry_1hd(tracer, active_species_idx, cp_dry, fact)
      use cam_abortutils,  only: endrun
      use string_utils,    only: int2str
      use physconst,       only: cpair
      use air_composition, only: dry_air_species_num

      ! Dummy arguments
      ! tracer: tracer array
      real(r8),           intent(in)  :: tracer(:,:,:)
      integer,            intent(in)  :: active_species_idx(:)
      ! fact: optional dry pressure level thickness
      real(r8), optional, intent(in)  :: fact(:,:)
      ! cp_dry: dry air heat capacity under constant pressure
      real(r8),           intent(out) :: cp_dry(:,:)

      ! Local variables
      integer  :: idx, jdx, kdx , m_cnst, qdx
      ! factor: dry pressure level thickness
      real(r8) :: factor(SIZE(cp_dry, 1), SIZE(cp_dry, 2))
      real(r8) :: residual(SIZE(cp_dry, 1), SIZE(cp_dry, 2))
      real(r8) :: mmr
      character(len=*), parameter :: subname = 'get_cp_dry_1hd: '

      if (dry_air_species_num == 0) then
         ! dry air heat capacity not species dependent
         cp_dry = cpair
      else
         ! dry air heat capacity is species dependent
         if (present(fact)) then
            if (SIZE(fact, 1) /= SIZE(factor, 1)) then
               call endrun(subname//"SIZE mismatch in dimension 1 "//         &
                    int2str(SIZE(fact, 1))//' /= '//int2str(SIZE(factor, 1)))
            end if
            if (SIZE(fact, 2) /= SIZE(factor, 2)) then
               call endrun(subname//"SIZE mismatch in dimension 2 "//         &
                    int2str(SIZE(fact, 2))//' /= '//int2str(SIZE(factor, 2)))
            end if
            factor = fact(:,:)
         else
            factor = 1.0_r8
         end if

         cp_dry = 0.0_r8
         residual = 1.0_r8
         do qdx = 1, dry_air_species_num
            m_cnst = active_species_idx(qdx)
            do kdx = 1, SIZE(cp_dry, 2)
               do idx = 1, SIZE(cp_dry, 1)
                  mmr = tracer(idx, kdx, m_cnst) * factor(idx, kdx)
                  cp_dry(idx, jdx, kdx) = cp_dry(idx, jdx, kdx) +             &
                       (thermodynamic_active_species_cp(qdx) * mmr)
                  residual(idx, kdx) = residual(idx, kdx) - mmr
               end do
            end do
         end do
         qdx = 0 ! N2
         do kdx = 1, SIZE(cp_dry, 2)
            do idx = 1, SIZE(cp_dry, 1)
               cp_dry(idx, kdx) = cp_dry(idx, kdx) +                          &
                    (thermodynamic_active_species_cp(qdx) * residual(idx, kdx))
            end do
         end do
      end if
   end subroutine get_cp_dry_1hd

   !===========================================================================

   subroutine get_cp_dry_2hd(tracer, active_species_idx, cp_dry, fact)
      ! Version of get_cp_dry for arrays that have a second horizontal index

      ! Dummy arguments
      ! tracer: tracer array
      real(r8),           intent(in)  :: tracer(:,:,:,:)
      integer,            intent(in)  :: active_species_idx(:)
      ! fact:        optional dry pressure level thickness
      real(r8), optional, intent(in)  :: fact(:,:,:)
      ! cp_dry: dry air heat capacity under constant pressure
      real(r8),           intent(out) :: cp_dry(:,:,:)

      ! Local variable
      integer  :: jdx

      do jdx = 1, SIZE(cp_dry, 2)
         if (present(fact)) then
            call get_cp_dry(tracer(:,jdx,:,:), active_species_idx,            &
                 cp_dry(:,jdx,:), fact=fact(:,jdx,:))
         else
            call get_cp_dry(tracer(:,jdx,:,:), active_species_idx,            &
                 cp_dry(:,jdx,:))
         end if
      end do

   end subroutine get_cp_dry_2hd

   !===========================================================================

   !
   !***************************************************************************
   !
   ! get_cp: Compute generalized heat capacity at constant pressure
   !
   !***************************************************************************
   !
   subroutine get_cp_1hd(tracer, inv_cp, cp, dp_dry, active_species_idx_dycore)
      use cam_abortutils,  only: endrun
      use string_utils,    only: int2str
      use air_composition, only: dry_air_species_num

      ! Dummy arguments
      ! tracedr: Tracer array
      real(r8),           intent(in)  :: tracer(:,:,:)
      real(r8), optional, intent(in)  :: dp_dry(:,:)
      ! inv_cp: output inverse cp instead of cp
      logical,            intent(in)  :: inv_cp
      real(r8),           intent(out) :: cp(:,:)
      ! active_species_idx_dycore: array of indicies for index of
      !    thermodynamic active species in dycore tracer array
      !    (if different from physics index)
      integer, optional,  intent(in)  :: active_species_idx_dycore(:)

      ! Local variables
      integer  :: qdx, itrac
      real(r8) :: sum_species(SIZE(cp, 1), SIZE(cp, 2))
      real(r8) :: sum_cp(SIZE(cp, 1), SIZE(cp, 2))
      real(r8) :: factor(SIZE(cp, 1), SIZE(cp, 2))
      integer  :: idx_local(thermodynamic_active_species_num)
      character(len=*), parameter :: subname = 'get_cp_1hd: '

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

      if (present(dp_dry)) then
         factor = 1.0_r8 / dp_dry
      else
         factor = 1.0_r8
      end if

      sum_species = 1.0_r8 ! all dry air species sum to 1
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
         itrac = idx_local(qdx)
         sum_species(:,:,:) = sum_species(:,:,:) +                            &
              (tracer(:,:,:,itrac) * factor(:,:,:))
      end do

      if (dry_air_species_num == 0) then
         sum_cp = thermodynamic_active_species_cp(0)
      else
         call get_cp_dry(tracer, idx_local, sum_cp, fact=factor)
      end if
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
         itrac = idx_local(qdx)
         sum_cp(:,:,:) = sum_cp(:,:,:) +                                      &
              (thermodynamic_active_species_cp(qdx) * tracer(:,:,:,itrac) *   &
              factor(:,:,:))
      end do
      if (inv_cp) then
         cp = sum_species / sum_cp
      else
         cp = sum_cp / sum_species
      end if

   end subroutine get_cp_1hd

   !===========================================================================

   subroutine get_cp_2hd(tracer, inv_cp, cp, dp_dry, active_species_idx_dycore)
      ! Version of get_cp for arrays that have a second horizontal index
      use cam_abortutils, only: endrun
      use string_utils,   only: int2str

      ! Dummy arguments
      ! tracedr: Tracer array
      real(r8),           intent(in)  :: tracer(:,:,:,:)
      real(r8), optional, intent(in)  :: dp_dry(:,:,:)
      ! inv_cp: output inverse cp instead of cp
      logical,            intent(in)  :: inv_cp
      real(r8),           intent(out) :: cp(:,:,:)
      ! active_species_idx_dycore: array of indicies for index of
      !    thermodynamic active species in dycore tracer array
      !    (if different from physics index)
      integer, optional,  intent(in)  :: active_species_idx_dycore(:)

      ! Local variables
      integer  :: jdx
      integer  :: idx_local(thermodynamic_active_species_num)
      character(len=*), parameter :: subname = 'get_cp_2hd: '

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
      do jdx = 1, SIZE(cp, 2)
         if (present(dp_dry)) then
            call get_cp(tracer(:, jdx, :, :), inv_cp, cp(:, jdx, :),          &
                 dp_dry=dp_dry(:, jdx, :), active_species_idx_dycore=idx_local)
         else
            call get_cp(tracer(:, jdx, :, :), inv_cp, cp(:, jdx, :),          &
                 active_species_idx_dycore=idx_local)
         end if
      end do

   end subroutine get_cp_2hd

   !===========================================================================

   !***************************************************************************
   !
   ! get_R_dry: Compute generalized dry air gas constant R
   !
   !***************************************************************************
   !
   subroutine get_R_dry_1hd(tracer, active_species_idx_dycore, R_dry, fact)
      use physconst,       only: rair
      use air_composition, only: dry_air_species_num

      ! tracer: tracer array
      real(r8),           intent(in)  :: tracer(:, :, :)
      ! active_species_idx_dycore: index of active species in tracer
      integer,            intent(in)  :: active_species_idx_dycore(:)
      ! R_dry: dry air R
      real(r8),           intent(out) :: R_dry(:, :)
      ! fact:   optional factor for converting tracer to dry mixing ratio
      real(r8), optional, intent(in)  :: fact(:, :)

      ! Local variables
      integer  :: idx, kdx, m_cnst, qdx
      real(r8) :: factor(SIZE(R_dry, 1), SIZE(R_dry, 2))
      real(r8) :: residual(SIZE(R_dry, 1), SIZE(R_dry, 2))
      real(r8) :: mmr

      if (dry_air_species_num == 0) then
         !
         ! dry air not species dependent
         !
         R_dry = rair
      else
         if (present(fact)) then
            factor = fact(:,:)
         else
            factor = 1.0_r8
         end if

         R_dry = 0.0_r8
         residual = 1.0_r8
         do qdx = 1, dry_air_species_num
            m_cnst = active_species_idx_dycore(qdx)
            do kdx = 1, SIZE(R_dry, 2)
               do idx = 1, SIZE(R_dry, 1)
                  mmr = tracer(idx, kdx, m_cnst) * factor(idx, kdx)
                  R_dry(idx, kdx) = R_dry(idx, kdx) +                         &
                       (thermodynamic_active_species_R(qdx) * mmr)
                  residual(idx, kdx) = residual(idx, kdx) - mmr
               end do
            end do
         end do
         !
         ! N2 derived from the others
         !
         qdx = 0
         do kdx = 1, SIZE(R_dry, 2)
            do idx = 1, SIZE(R_dry, 1)
               R_dry(idx, kdx) = R_dry(idx, kdx) +                            &
                    (thermodynamic_active_species_R(qdx) * residual(idx, kdx))
            end do
         end do
      end if
   end subroutine get_R_dry_1hd

   !===========================================================================

   subroutine get_R_dry_2hd(tracer, active_species_idx_dycore, R_dry, fact)
      ! Version of get_R_dry for arrays that have a second horizontal index

      ! tracer: tracer array
      real(r8),           intent(in)  :: tracer(:, :, :, :)
      ! active_species_idx_dycore: index of active species in tracer
      integer,            intent(in)  :: active_species_idx_dycore(:)
      ! R_dry: dry air R
      real(r8),           intent(out) :: R_dry(:, :, :)
      ! fact:   optional factor for converting tracer to dry mixing ratio
      real(r8), optional, intent(in)  :: fact(:, :, :)

      ! Local variable
      integer  :: jdx

      do jdx = 1, SIZE(R_dry, 2)
         if (present(fact)) then
            call get_R_dry(tracer(:, jdx, :), active_species_idx_dycore,      &
                 R_dry(:, jdx, :), fact=fact(:, jdx, :))
         else
            call get_R_dry(tracer(:, jdx, :), active_species_idx_dycore,      &
                 R_dry(:, jdx, :))
         end if
      end do

   end subroutine get_R_dry_2hd

   !
   !***************************************************************************
   !
   ! get_R: Compute generalized R
   !
   !***************************************************************************
   !
   subroutine get_R_1hd(tracer, active_species_idx, R, fact)
      use cam_abortutils,  only: endrun
      use string_utils,    only: int2str
      use air_composition, only: dry_air_species_num

      ! Dummy arguments
      ! tracer: !tracer array
      real(r8), intent(in)  :: tracer(:, :, :)
      ! active_species_idx: index of active species in tracer
      integer,  intent(in)  :: active_species_idx(:)
      ! R: generalized gas constant
      real(r8), intent(out) :: R(:, :)
      ! fact: optional factor for converting tracer to dry mixing ratio
      real(r8), optional, intent(in) :: fact(:, :)

      ! Local variables
      integer  :: qdx, itrac
      real(r8) :: factor(SIZE(R, 1), SIZE(R, 2))
      real(r8) :: sum_species(SIZE(R, 1), SIZE(R, 2))
      integer  :: idx_local(thermodynamic_active_species_num)

      character(len=*), parameter :: subname = 'get_R_1hd: '

      if (present(fact)) then
         if (SIZE(fact, 1) /= SIZE(factor, 1)) then
            call endrun(subname//"SIZE mismatch in dimension 1 "//            &
                 int2str(SIZE(fact, 1))//' /= '//int2str(SIZE(factor, 1)))
         end if
         if (SIZE(fact, 2) /= SIZE(factor, 2)) then
            call endrun(subname//"SIZE mismatch in dimension 2 "//            &
                 int2str(SIZE(fact, 2))//' /= '//int2str(SIZE(factor, 2)))
         end if
         call get_R_dry(tracer, active_species_idx, R, fact=fact)
         factor = fact(:,:,:)
      else
         call get_R_dry(tracer, active_species_idx, R)
         factor = 1.0_r8
      end if
      idx_local = active_species_idx
      sum_species = 1.0_r8 ! all dry air species sum to 1
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
         itrac = idx_local(qdx)
         sum_species(:,:,:) = sum_species(:,:,:) +                            &
              (tracer(:,:,:,itrac) * factor(:,:,:))
      end do
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
         itrac = idx_local(qdx)
         R(:,:,:) = R(:,:,:) +                                                &
              (thermodynamic_active_species_R(qdx) * tracer(:,:,:,itrac) *    &
              factor(:,:,:))
      end do
      R = R / sum_species
   end subroutine get_R_1hd

   !===========================================================================

   subroutine get_R_2hd(tracer, active_species_idx, R, fact)

      ! Dummy arguments
      ! tracer: !tracer array
      real(r8), intent(in)  :: tracer(:, :, :, :)
      ! active_species_idx: index of active species in tracer
      integer,  intent(in)  :: active_species_idx(:)
      ! R: generalized gas constant
      real(r8), intent(out) :: R(:, :, :)
      ! fact: optional factor for converting tracer to dry mixing ratio
      real(r8), optional, intent(in) :: fact(:, :, :)

      ! Local variable
      integer  :: jdx

      do jdx = 1, SIZE(R, 2)
         if (present(fact)) then
            call get_R(tracer(:, jdx, :), active_species_idx_dycore,          &
                 R(:, jdx, :), fact=fact(:, jdx, :))
         else
            call get_R(tracer(:, jdx, :), active_species_idx_dycore,          &
                 R(:, jdx, :))
         end if
      end do

   end subroutine get_R_2hd

   !===========================================================================

   subroutine get_thermal_energy_1hd(tracer_mass, temp, dp_dry,               &
        thermal_energy, active_species_idx_dycore)
      use air_composition, only: dry_air_species_num

      !
      !***********************************************************************
      !
      ! Compute thermal energy = cp*T*dp, where dp is pressure level thickness,
      !    cp is generalized cp and T temperature
      !
      ! Note: tracer is in units of m*dp_dry ("mass")
      !
      !***********************************************************************
      !
      ! Dummy arguments
      ! tracer_mass: tracer array (mass weighted)
      real(r8),          intent(in)  :: tracer_mass(:,:,:,:)
      ! temp: temperature
      real(r8),          intent(in)  :: temp(:,:,:)
      ! dp_dry: dry presure level thickness
      real(r8),          intent(in)  :: dp_dry(:,:,:)
      ! thermal_energy: thermal energy in each column: sum cp*T*dp
      real(r8),          intent(out) :: thermal_energy(:,:,:)
      !
      ! active_species_idx_dycore:
      !    array of indicies for index of thermodynamic active species in
      !    dycore tracer array (if different from physics index)
      !
      integer, optional, intent(in)  :: active_species_idx_dycore(:)

      ! Local vars
      integer                     :: qdx, itrac
      character(len=*), parameter :: subname = 'get_thermal_energy: '

      !
      ! "mass-weighted" cp (dp must be dry)
      !
      if (dry_air_species_num == 0) then
         thermal_energy(:,:,:) = thermodynamic_active_species_cp(0) *         &
              dp_dry(:,:,:)
      else
         if (present(active_species_idx_dycore)) then
            call get_cp_dry(tracer_mass, active_species_idx_dycore,           &
                 thermal_energy, fact=1.0_r8/dp_dry(:,:,:))
         else
            call get_cp_dry(tracer_mass, thermodynamic_active_species_idx,    &
                 thermal_energy, fact=1.0_r8/dp_dry(:,:,:))
         end if
         thermal_energy(:,:,:) = thermal_energy(:,:,:) * dp_dry(:,:,:)
      end if
      !
      ! tracer is in units of m*dp ("mass"), where:
      !    m is the dry mixing ratio
      !    dp is the dry pressure level thickness
      !
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
         if (present(active_species_idx_dycore)) then
            itrac = active_species_idx_dycore(qdx)
         else
            itrac = thermodynamic_active_species_idx(qdx)
         end if
         thermal_energy(:,:,:) = thermal_energy(:,:,:) +                      &
              (thermodynamic_active_species_cp(qdx) * tracer_mass(:,:,:,itrac))
      end do
      thermal_energy(:,:,:) = thermal_energy(:,:,:) * temp(:,:,:)

   end subroutine get_thermal_energy_1hd

   !===========================================================================

   subroutine get_thermal_energy_2hd(tracer_mass, temp, dp_dry,               &
        thermal_energy, active_species_idx_dycore)
      !
      !***********************************************************************
      !
      ! Compute thermal energy = cp*T*dp, where dp is pressure level thickness,
      !    cp is generalized cp and T temperature
      !
      ! Note: tracer is in units of m*dp_dry ("mass")
      !
      !***********************************************************************
      !
      ! Dummy arguments
      ! tracer_mass: tracer array (mass weighted)
      real(r8),          intent(in)  :: tracer_mass(:,:,:,:)
      ! temp: temperature
      real(r8),          intent(in)  :: temp(:,:,:)
      ! dp_dry: dry presure level thickness
      real(r8),          intent(in)  :: dp_dry(:,:,:)
      ! thermal_energy: thermal energy in each column: sum cp*T*dp
      real(r8),          intent(out) :: thermal_energy(:,:,:)
      !
      ! active_species_idx_dycore:
      !    array of indicies for index of thermodynamic active species in
      !    dycore tracer array (if different from physics index)
      !
      integer, optional, intent(in)  :: active_species_idx_dycore(:)

      ! Local variables
      integer                     :: jdx
      character(len=*), parameter :: subname = 'get_thermal_energy_2hd: '

      do jdx = 1, SIZE(tracer_mass, 2)
         call get_thermal_energy(tracer_mass(:, jdx, :), temp(:, jdx, :),     &
              dp_dry(:, jdx, :), thermal_energy(:, jdx, :),                   &
              active_species_idx_dycore=active_species_idx_dycore)
      end do

   end subroutine get_thermal_energy_2hd

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
      use air_composition, only: dry_air_species_num

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
      real(r8)                    :: sum_species
      real(r8)                    :: factor
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

      if (present(dp_dry)) then
         factor = 1.0_r8 / dp_dry
      else
         factor = 1.0_r8
      end if

      sum_species = 1.0_r8 ! All dry air species sum to 1
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
         itrac = idx_local(qdx)
         sum_species(:, :) = sum_species(:, :) +                              &
              (tracer(:, :, itrac) * factor(:, :))
      end do

      call get_R_dry(tracer, idx_local, Rd, fact=factor)
      t_v(:, :) = Rd(:, :)
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
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
        sum_species,dp_dry)
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

      ! Local variables
      real(r8) :: factor(SIZE(tracer, 1), SIZE(tracer, 2))
      integer  :: qdx, itrac

      if (present(dp_dry)) then
         factor = 1.0_r8 / dp_dry(:,:)
      else
         factor = 1.0_r8
      end if
      sum_species = 1.0_r8 ! all dry air species sum to 1
      do qdx = dry_air_species_num+1, thermodynamic_active_species_num
         itrac = active_species_idx(qdx)
         sum_species(:,:) = sum_species(:,:) + (tracer(:,:,itrac) * factor(:,:))
      end do
   end subroutine get_sum_species_1hd

   !===========================================================================

   subroutine get_sum_species_2hd(tracer, active_species_idx,                 &
        sum_species,dp_dry)

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

      ! Local variable
      integer                     :: jdx

      do jdx = 1, SIZE(tracer, 2)
         if (present(dp_dry)) then
            call get_sum_species(tracer(:, jdx, :, :), active_species_idx,    &
                 dp_dry=dp_dry(:, jdx, :), sum_species(:, jdx, :))
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

     integer :: i, j, k, m_cnst, nq

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

   !===========================================================================

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

   !===========================================================================

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

     if (dycore_is('LR').or.dycore_is('FV3')) then
       do k=k0,k1
         pmid(:,:,k) = dp(:,:,k) / (log(pint_local(:,:,k+1))-log(pint_local(:,:,k)))
       end do
     else
       do k=k0,k1
         pmid(:,:,k) = 0.5_r8*(pint_local(:,:,k)+pint_local(:,:,k+1))
       end do
     end if
     if (present(pint)) pint=pint_local
   end subroutine get_pmid_from_dp

   !===========================================================================

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
     ! compute kappa = Rd / cpd
     !
     if (mixing_ratio==1) then
       call get_kappa_dry(i0,i1,j0,j1,1,nlev,nlev,ntrac,tracer,active_species_idx,kappa_dry)
     else
       call get_kappa_dry(i0,i1,j0,j1,1,nlev,nlev,ntrac,tracer,active_species_idx,kappa_dry,1.0_r8 / dp_dry)
     end if
     if (inv_exner) then
       exner(:,:,:) = (p00 / pmid(:,:,:))**kappa_dry(:,:,:)
     else
       exner(:,:,:) = (pmid(:,:,:) / p00)**kappa_dry(:,:,:)
     end if
     if (present(poverp0)) poverp0=pmid(:,:,:) / p00
   end subroutine get_exner

   !===========================================================================

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

   !===========================================================================

   !****************************************************************************************************************
   !
   ! Compute geopotential from dry pressure level thichkness, water tracers, model top pressure and temperature
   !
   !****************************************************************************************************************
   !
   subroutine get_gz_from_dp_dry_ptop_temp(i0,i1,j0,j1,nlev,ntrac,tracer,mixing_ratio,active_species_idx, &
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
                      R_dry,fact=1.0_r8 / dp_dry)
     end if
     call get_gz_given_dp_Tv_Rdry(i0,i1,j0,j1,nlev,dp_local,T_v_local,R_dry,phis,ptop,gz,pmid_local)

     if (present(pmid)) pmid=pmid_local
     if (present(T_v))  T_v=T_v_local
     if (present(dp))   dp=dp_local
  end subroutine get_gz_from_dp_dry_ptop_temp

   !===========================================================================

   !***************************************************************************
   !
   ! Compute geopotential from pressure level thickness and virtual temperature
   !
   !***************************************************************************
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
     if (dycore_is('LR').or.dycore_is('FV3')) then
       do k=nlev,1,-1
         Rdry_tv(:,:) = R_dry(:,:,k)*T_v(:,:,k)
         gz(:,:,k) = gzh(:,:)+Rdry_tv(:,:)*(1.0_r8-pint(:,:,k) / pmid_local(:,:,k))
         gzh(:,:)  = gzh(:,:) + Rdry_tv(:,:)*(log(pint(:,:,k+1))-log(pint(:,:,k)))
       end do
     else
       do k=nlev,1,-1
         Rdry_tv(:,:) = R_dry(:,:,k)*T_v(:,:,k)
         gz(:,:,k) = gzh(:,:)+Rdry_tv(:,:)*0.5_r8*dp(:,:,k) / pmid_local(:,:,k)
         gzh(:,:)  = gzh(:,:) + Rdry_tv(:,:)*dp(:,:,k) / pmid_local(:,:,k)
       end do
     end if
     if (present(pmid)) pmid=pmid_local
   end subroutine get_gz_given_dp_Tv_Rdry

   !===========================================================================

   !***************************************************************************
   !
   ! Compute Richardson number at cell interfaces (half levels)
   !
   !***************************************************************************
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
       Richardson_number(:,:,k) = (gz(:,:,km1)-gz(:,:,k))*(pt1-pt2) / ( 0.5_r8*(pt1+pt2)*        &
            ((v(:,:,1,km1)-v(:,:,1,k))**2+(v(:,:,2,km1)-v(:,:,2,k))**2+ustar2) )
     end do
   end subroutine get_Richardson_number

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

   !===========================================================================

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
         factor(i0:i1,j0:j1,k0:k1) = fact(i0:i1,j0:j1,k0:k1)
       else
         factor(i0:i1,j0:j1,k0:k1) = 1.0_r8
       endif

       mbarv = 0.0_r8
       residual = 1.0_r8
       do nq=1,dry_air_species_num
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
       nq = 0 ! N2
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

   !===========================================================================

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

   !===========================================================================

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

   !===========================================================================

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

   !===========================================================================

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
               do icnst=1,dry_air_species_num
                 ispecies = idx_local(icnst)
                 mm       = 0.5_r8*(tracer(i,j,k,ispecies)*factor(i,j,k)+tracer(i,j,k-1,ispecies)*factor(i,j,k-1))
                 kmvis(i,j,k) = kmvis(i,j,k)+thermodynamic_active_species_kv(icnst)* &
                                             thermodynamic_active_species_mwi(icnst)*mm
                 kmcnd(i,j,k) = kmcnd(i,j,k)+thermodynamic_active_species_kc(icnst)* &
                                             thermodynamic_active_species_mwi(icnst)*mm
                 residual         = residual - mm
               end do
               icnst=0 ! N2
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
         do k=1,k1
           do j=j0,j1
             do i=i0,i1
               kmvis(i,j,k) = 0.0_r8
               kmcnd(i,j,k) = 0.0_r8
               residual = 1.0_r8
               do icnst=1,dry_air_species_num-1
                 ispecies = idx_local(icnst)
                 mm       = tracer(i,j,k,ispecies)*factor(i,j,k)
                 kmvis(i,j,k) = kmvis(i,j,k)+thermodynamic_active_species_kv(icnst)* &
                                thermodynamic_active_species_mwi(icnst)*mm
                 kmcnd(i,j,k) = kmcnd(i,j,k)+thermodynamic_active_species_kc(icnst)* &
                                thermodynamic_active_species_mwi(icnst)*mm
                 residual     = residual - mm
               end do
               icnst=dry_air_species_num
               kmvis(i,j,k) = kmvis(i,j,k)+thermodynamic_active_species_kv(icnst)* &
                              thermodynamic_active_species_mwi(icnst)*residual
               kmcnd(i,j,k) = kmcnd(i,j,k)+thermodynamic_active_species_kc(icnst)* &
                              thermodynamic_active_species_mwi(icnst)*residual

               kmvis(i,j,k) = kmvis(i,j,k)*mbarv(i,j,k)*temp(i,j,k)**kv4*1.e-7_r8
               kmcnd(i,j,k) = kmcnd(i,j,k)*mbarv(i,j,k)*temp(i,j,k)**kc4*1.e-5_r8
             enddo
           enddo
         end do
       else
         call endrun('get_molecular_diff_coef: get_at_interfaces must be 0 or 1')
       end if
     end if
   end subroutine get_molecular_diff_coef

   !===========================================================================

   !***************************************************************************
   !
   ! compute reference vertical profile of density, molecular diffusion and thermal conductivity
   !
   !***************************************************************************
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

   !==========================================================================

   !
   !***************************************************************************
   !
   ! cam_thermo_calc_kappav: update species dependent kappa for FV dycore
   !
   !***************************************************************************
   !
   subroutine cam_thermo_calc_kappav_2hd(tracer, kappav, cpv)
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
      !$omp parallel do private(i,j,k)
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
   subroutine get_hydrostatic_energy_1hd(tracer, pdel, cp_or_cv, U, V, T,     &
        vcoord, ps, phis, z_mid, dycore_idx, qidx, te, se, ke,                &
        wv, H2O, liq, ice)

      use cam_logfile,     only: iulog
      use dyn_tests_utils, only: vc_height, vc_moist_pressure, vc_dry_pressure
      use air_composition, only: wv_idx

      ! Dummy arguments
      ! tracer: tracer mixing ratio
      real(r8), intent(in)            :: tracer(:,:,:)
      ! pdel: pressure level thickness
      real(r8), intent(in)            :: pdel(:,:)
      ! cp_or_cv: dry air heat capacity under constant pressure or
      !           constant volume (depends on vcoord)
      real(r8), intent(in)            :: cp_or_cv(:,:)
      real(r8), intent(in)            :: U(:,:)
      real(r8), intent(in)            :: V(:,:)
      real(r8), intent(in)            :: T(:,:)
      integer,  intent(in)            :: vcoord ! vertical coordinate
      real(r8), intent(in),  optional :: ps(:)
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
      ! SE: vertically integrated internal+geopotential energy
      real(r8), intent(out), optional :: se (:)
      ! WV: vertically integrated water vapor
      real(r8), intent(out), optional :: wv (:)
      ! liq: vertically integrated liquid
      real(r8), intent(out), optional :: liq(:)
      ! ice: vertically integrated ice
      real(r8), intent(out), optional :: ice(:)

      ! Local variables
      real(r8) :: ke_vint(SIZE(tracer, 1))  ! Vertical integral of KE
      real(r8) :: se_vint(SIZE(tracer, 1))  ! Vertical integral of SE
      real(r8) :: wv_vint(SIZE(tracer, 1))  ! Vertical integral of wv
      real(r8) :: liq_vint(SIZE(tracer, 1)) ! Vertical integral of liq
      real(r8) :: ice_vint(SIZE(tracer, 1)) ! Vertical integral of ice
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

      select case (vcoord)
      case(vc_moist_pressure, vc_dry_pressure)
         if ((.not. present(ps)) .or. (.not. present(phis))) then
            write(iulog, *) subname, ' ps and phis must be present for ',     &
                 'moist/dry pressure vertical coordinate'
            call endrun(subname//':  ps and phis must be present for '//      &
                 'moist/dry pressure vertical coordinate')
         end if
         ke_vint = 0._r8
         se_vint = 0._r8
         wv_vint = 0._r8
         do kdx = 1, SIZE(tracer, 2)
            do idx = 1, SIZE(tracer, 1)
               ke_vint(idx) = ke_vint(idx) + (pdel(idx, kdx) *                &
                    0.5_r8 * (U(idx, kdx)**2 + V(idx, kdx)**2) / gravit)
               se_vint(i,j) = se_vint(i,j) + (T(idx, kdx) *                   &
                    cp_or_cv(idx, kdx) * pdel(idx, kdx) / gravit)
               wv_vint(idx) = wv_vint(idx) + (tracer(idx, kdx, wvidx) *       &
                    pdel(idx, kdx) / gravit)
            end do
         end do
         do idx = 1, SIZE(tracer, 1)
            se_vint(idx) = se_vint(idx) + (phis(idx) * ps(idx) / gravit)
         end do
      case(vc_height)
         if (.not. present(z_mid)) then
            write(iulog, *) subname,                                          &
                 ' z_mid must be present for height vertical coordinate'
            call endrun(subname//': z_mid must be present for height '//      &
                 'vertical coordinate')
         end if
         ke_vint = 0._r8
         se_vint = 0._r8
         wv_vint = 0._r8
         do kdx = 1, SIZE(tracer, 2)
            do i = 1, SIZE(tracer, 1)
               ke_vint(idx) = ke_vint(idx) + (pdel(idx, kdx) *                &
                    0.5_r8 * (U(idx, kdx)**2 + V(idx, kdx)**2) / gravit)
               se_vint(idx) = se_vint(idx) + (T(idx, kdx) *                   &
                    cp_or_cv(idx, kdx) * pdel(idx, kdx) / gravit)
               ! z_mid is height above ground
               se_vint(idx) = se_vint(idx) + (z_mid(idx, kdx) +               &
                    (phis(idx) / gravit) * pdel(idx, kdx))
               wv_vint(idx) = wv_vint(idx) + (tracer(idx, kdx, wvidx) *       &
                    pdel(idx, kdx) / gravit)
            end do
         end do
      case default
         write(iulog, *) subname, ' vertical coordinate not supported: ', vcoord
         call endrun(subname//': vertical coordinate not supported')
      end select
      if (present(te)) then
         te  = se_vint + ke_vint
      end if
      if (present(se)) then
         se = se_vint
      end if
      if (present(ke)) then
         ke = ke_vint
      end if
      if (present(wv)) then
         wv = wv_vint
      end if
      !
      ! vertical integral of total liquid water
      !
      liq_vint = 0._r8
      do qdx = 1, thermodynamic_active_species_liq_num
         do kdx = 1, SIZE(tracer, 2)
            do idx = 1, SIZE(tracer, 1)
               liq_vint(qdx) = liq_vint(qdx) + (pdel(idx, kdx) *              &
                    tracer(idx, kdx, species_liq_idx(qdx)) / gravit)
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
               ice_vint(qdx) = ice_vint(qdx) + (pdel(idx, kdx) *              &
                    tracer(idx, kdx, species_ice_idx(qdx))  / gravit)
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

   !===========================================================================

end module cam_thermo
