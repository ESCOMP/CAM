      module chem_mods
!--------------------------------------------------------------
! ... Basic chemistry parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod,    only : r8 => shr_kind_r8, shr_kind_cl
      use constituents,    only : pcnst
      implicit none
      save

      INTEGER, PARAMETER   :: nTracersMax = 267    ! Must be equal to chem_nadv
      INTEGER              :: nTracers
      REAL(r8)             :: ref_MMR(pcnst)

      CHARACTER(LEN=shr_kind_cl) :: tracerNames(nTracersMax)
      CHARACTER(LEN=shr_kind_cl) :: tracerLongNames(nTracersMax)

      ! Index of first constituent
      INTEGER              :: iFirstCnst

      ! Short-lived species (i.e. not advected)
      INTEGER, PARAMETER   :: nSlsMax = 500        ! UNadvected species only
      INTEGER              :: nSls

      CHARACTER(LEN=shr_kind_cl) :: slsNames(nSlsMax)
      CHARACTER(LEN=shr_kind_cl) :: slsLongnames(nSlsMax)

      ! Mapping between constituents and GEOS-Chem tracers
      INTEGER              :: map2GC(pcnst)
      INTEGER              :: map2GCinv(nTracersMax)
      INTEGER              :: map2GC_Sls(nSlsMax)

      ! Mapping constituent onto chemical species (as listed in solsym)
      INTEGER              :: mapCnst(pcnst)

      ! Aerosols
      INTEGER, PARAMETER   :: nAerMax = 35
      INTEGER              :: nAer
      CHARACTER(LEN=16)    :: aerNames(nAerMax)
      REAL(r8)             :: aerAdvMass(nAerMax)

      !-----------------------------
      ! Aerosol index mapping
      !-----------------------------
      ! map2MAM4 maps aerNames onto the GEOS-Chem Species array such
      ! that
      ! State_Chm%Species(1,:,:,map2MAM4(:,:)) = state%q(:,:,MAM4_Indices)
      INTEGER, ALLOCATABLE :: map2MAM4(:,:)

      !-----------------------------
      ! Dry deposition index mapping
      !-----------------------------
      ! drySpc_ndx maps drydep_list onto tracerNames such that
      ! tracerNames(drySpc_ndx(:)) = drydep_list(:)
      INTEGER, ALLOCATABLE :: drySpc_ndx(:)

      ! map2GC_dryDep maps drydep_list onto the GEOS-Chem dry deposition
      ! velocity arrays such that
      ! State_Chm%DryDepVel(1,:,map2GC_dryDep(:)) = cam_in%depVel(:,:)
      INTEGER, ALLOCATABLE :: map2GC_dryDep(:)

      INTEGER, PARAMETER :: phtcnt = 40, & ! number of photolysis reactions
                            rxntot = 212, & ! number of total reactions
                            gascnt = 172, & ! number of gas phase reactions
                            nabscol = 2, & ! number of absorbing column densities
                            gas_pcnst = 269, & ! number of "gas phase" species (same as solsym length)
                                               ! Includes GC advected species (233), MAM aerosols (33),
                                               ! and CO2 (1), as well as any non-advected species added
                                               ! to solsym and mo_sim_dat.F90.
                            nfs = 6, & ! number of "fixed" species
                            relcnt = 0, & ! number of relationship species
                            grpcnt = 0, & ! number of group members
                            nzcnt = 824, & ! number of non-zero matrix entries
                            extcnt = 34, & ! number of species with external forcing, aka 3-D emissions
                            clscnt1 = 8, & ! number of species in explicit class
                            clscnt2 = 0, & ! number of species in hov class
                            clscnt3 = 0, & ! number of species in ebi class
                            clscnt4 = 95, & ! number of species in implicit class
                            clscnt5 = 0, & ! number of species in rodas class
                            indexm = 1, & ! index of total atm density in invariant array
                            indexh2o = 4, & ! index of water vapor density
                            clsze = 1, & ! loop length for implicit chemistry
                            rxt_tag_cnt = 0, & ! number of tagged reactions (unused in GEOS-Chem)
                            enthalpy_cnt = 0, &
                            nslvd = 86  ! number of short-lived (non-advected) species
      integer :: clscnt(5) = 0
      integer :: cls_rxt_cnt(4,5) = 0
      integer :: clsmap(gas_pcnst,5) = 0
      integer :: permute(gas_pcnst,5) = 0
      integer :: diag_map(clscnt4) = 0
      real(r8) :: adv_mass(gas_pcnst) = 0._r8
      real(r8) :: crb_mass(gas_pcnst) = 0._r8
      real(r8) :: fix_mass(max(1,nfs))
      real(r8), allocatable :: cph_enthalpy(:)
      integer, allocatable :: cph_rid(:)
      integer, allocatable :: num_rnts(:)
      integer, allocatable :: rxt_tag_map(:)
      real(r8), allocatable :: pht_alias_mult(:,:)
      character(len=16), allocatable :: rxt_tag_lst(:)
      character(len=16), allocatable :: pht_alias_lst(:,:)
      character(len=16) :: inv_lst(max(1,nfs))
      character(len=16) :: extfrc_lst(max(1,extcnt))
      logical :: frc_from_dataset(max(1,extcnt))
      logical :: is_vector
      logical :: is_scalar
      character(len=shr_kind_cl), allocatable :: slvd_lst(:)

      ! Mapping between chemical species and GEOS-Chem species/other tracers
      INTEGER              :: map2chm(gas_pcnst)

      end module chem_mods
