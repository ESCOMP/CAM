      module chem_mods
!--------------------------------------------------------------
! ... Basic chemistry parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod,    only : r8 => shr_kind_r8
      use constituents,    only : pcnst
      implicit none
      save

      INTEGER, PARAMETER   :: nTracersMax = 240    ! Must be equal to chem_nadv
      INTEGER              :: nTracers
      CHARACTER(LEN=255)   :: tracerNames(nTracersMax)
      CHARACTER(LEN=255)   :: tracerLongNames(nTracersMax)
      REAL(r8)             :: adv_Mass(nTracersMax)
      REAL(r8)             :: MWRatio(nTracersMax)
      REAL(r8)             :: ref_MMR(nTracersMax)

      ! Short-lived species (i.e. not advected)
      INTEGER, PARAMETER   :: nSlsMax = 500        ! UNadvected species only
      INTEGER              :: nSls    
      CHARACTER(LEN=255)   :: slsNames(nSlsMax)
      CHARACTER(LEN=255)   :: slsLongnames(nSlsMax)
      REAL(r8)             :: sls_Ref_MMR(nSlsMax)
      REAL(r8)             :: slsMWRatio(nSlsMax)

      ! Mapping between constituents and GEOS-Chem tracers
      INTEGER              :: map2GC(pcnst)
      INTEGER              :: map2GC_Sls(nSlsMax)

      ! Mapping from constituents to raw index
      INTEGER              :: map2Idx(pcnst)

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
                            gas_pcnst = 103, & ! number of "gas phase" species
                            nfs = 6, & ! number of "fixed" species
                            relcnt = 0, & ! number of relationship species
                            grpcnt = 0, & ! number of group members
                            nzcnt = 824, & ! number of non-zero matrix entries
                            extcnt = 4, & ! number of species with external forcing
                            clscnt1 = 8, & ! number of species in explicit class
                            clscnt2 = 0, & ! number of species in hov class
                            clscnt3 = 0, & ! number of species in ebi class
                            clscnt4 = 95, & ! number of species in implicit class
                            clscnt5 = 0, & ! number of species in rodas class
                            indexm = 1, & ! index of total atm density in invariant array
                            indexh2o = 4, & ! index of water vapor density
                            clsze = 1, & ! loop length for implicit chemistry
                            rxt_tag_cnt = 95, &
                            enthalpy_cnt = 0
!                            nslvd = 0
      integer :: clscnt(5) = 0
      integer :: cls_rxt_cnt(4,5) = 0
      integer :: clsmap(gas_pcnst,5) = 0
      integer :: permute(gas_pcnst,5) = 0
      integer :: diag_map(clscnt4) = 0
      !real(r8) :: adv_mass(gas_pcnst) = 0._r8
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
!      character(len=16) :: slvd_lst(max(1,nslvd))
      integer :: nslvd
      character(len=255), allocatable :: slvd_lst(:)
      real(r8), allocatable :: slvd_ref_mmr(:)
      end module chem_mods
