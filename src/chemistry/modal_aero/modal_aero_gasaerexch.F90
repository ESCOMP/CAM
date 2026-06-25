! Portable code for modal aerosol gas-aerosol exchange.
! RCE 07.04.13:  Adapted from MIRAGE2 code
module modal_aero_gasaerexch
  use shr_kind_mod,  only: r8 => shr_kind_r8

  implicit none
  private
  save

  public :: modal_aero_gasaerexch_init
  public :: modal_aero_gasaerexch_run

  ! Primary-carbon aging (pcage) configuration: species are transferred from
  ! the primary-carbon mode (modefrm) to the accumulation mode (modetoo) when
  ! enough sulfate monolayers coat the particle.
  integer, protected, public :: maxspec_pcage     ! max number of species that can be aged

  integer, protected, public :: modefrm_pcage     ! source mode index for aging transfer
  integer, protected, public :: nspecfrm_pcage    ! number of species transferred during aging

  integer, protected, allocatable, public :: lspecfrm_pcage(:) ! pcnst indices of species in source mode
  integer, protected, allocatable, public :: lspectoo_pcage(:) ! pcnst indices of corresponding species in dest mode

  real(r8), parameter, public :: n_so4_monolayers_pcage = 8.0_r8

  ! number of so4(+nh4) monolayers needed to "age" a carbon particle
  ! thickness of the so4 monolayers (m)
  ! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3,
  !    --> 1 mol so4(+nh4)  = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
  ! aging criterion is approximate so do not try to distinguish
  !    sulfuric acid, bisulfate, ammonium sulfate
  real(r8), parameter, public :: &
              dr_so4_monolayers_pcage = n_so4_monolayers_pcage * 4.76e-10_r8

  ! this factor converts an soa volume to a volume of so4(+nh4)
  ! having same hygroscopicity as the soa
  real(r8), protected, allocatable, public :: soa_equivso4_factor(:)

  ! Private module-level storage:

  ! Mode configuration
  integer :: ntot_amode_m, nsoa_m, npoa_m, nspec_max_m
  integer, allocatable :: nspec_amode_m(:)

  ! Species indices in pcnst-space (set by _init, converted to vmr-space in _run)
  integer :: idx_h2so4_m, idx_nh3_m, idx_msa_m
  integer, allocatable :: idx_soag_m(:)
  integer, allocatable :: idx_so4_a_m(:), idx_nh4_a_m(:)
  integer, allocatable :: idx_soa_a_m(:,:), idx_pom_a_m(:,:)
  integer, allocatable :: idx_num_m(:), idx_mass_m(:,:)

  ! Mode metadata
  real(r8), allocatable :: alnsg_amode_m(:), sigmag_amode_m(:)
  real(r8), allocatable :: specmw_amode_m(:,:), specdens_amode_m(:,:)

  ! Flags
  logical :: do_nh4g_m, do_msag_m, do_soag_any_m
  logical, allocatable :: do_soag_m(:)

  ! Species presence in modes
  integer, allocatable :: ido_so4a_m(:), ido_nh4a_m(:), ido_soaa_m(:,:)
  integer :: ntot_soamode_m

  ! pcage and pcarbon
  integer :: modetoo_pcage
  integer :: modeptr_pcarbon_m

  ! Mass-to-volume conversion factors
  real (r8) :: fac_m2v_nh4, fac_m2v_so4
  real (r8), allocatable :: fac_m2v_soa(:)

  real (r8), allocatable :: fac_m2v_pcarbon(:)

  ! SOA/POA molecular weights from host model
  real(r8), allocatable :: mw_soa_host_m(:), mw_poa_host_m(:)

  ! Host-provided physical constants:
  real(r8) :: rair_m, mwdry_m
  real(r8) :: rgas_m

contains

subroutine modal_aero_gasaerexch_init( &
    ntot_amode, nsoa, npoa, nspec_max, &
    nspec_amode, &
    modeptr_pcarbon, modeptr_accum, &
    alnsg_amode, sigmag_amode, &
    specmw_amode, specdens_amode, spechygro, &
    idx_h2so4, idx_nh3, idx_msa, &
    idx_soag, &
    idx_so4_a, idx_nh4_a, &
    idx_soa_a, idx_pom_a, &
    idx_num, idx_mass, pcnst_in, &
    nspecfrm_pcage_in, &
    lspecfrm_pcage_in, lspectoo_pcage_in, &
    mw_soa_host, mw_poa_host, &
    rair, mwdry, r_universal, &
    errmsg, errflg)

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   !    initialize gas-aerosol exchange module
   !    store species indices and mode metadata
   !    compute aging/MW conversion factors
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   ! arguments
   integer,  intent(in) :: ntot_amode
   integer,  intent(in) :: nsoa
   integer,  intent(in) :: npoa
   integer,  intent(in) :: nspec_max
   integer,  intent(in) :: nspec_amode(:)
   integer,  intent(in) :: modeptr_pcarbon
   integer,  intent(in) :: modeptr_accum
   real(r8), intent(in) :: alnsg_amode(:)
   real(r8), intent(in) :: sigmag_amode(:)
   real(r8), intent(in) :: specmw_amode(:,:)
   real(r8), intent(in) :: specdens_amode(:,:)
   real(r8), intent(in) :: spechygro(:,:)
   integer,  intent(in) :: idx_h2so4
   integer,  intent(in) :: idx_nh3
   integer,  intent(in) :: idx_msa
   integer,  intent(in) :: idx_soag(:)
   integer,  intent(in) :: idx_so4_a(:)
   integer,  intent(in) :: idx_nh4_a(:)
   integer,  intent(in) :: idx_soa_a(:,:)
   integer,  intent(in) :: idx_pom_a(:,:)
   integer,  intent(in) :: idx_num(:)
   integer,  intent(in) :: idx_mass(:,:)
   integer,  intent(in) :: pcnst_in           ! total number of constituents (for range checks)
   integer,  intent(in) :: nspecfrm_pcage_in
   integer,  intent(in) :: lspecfrm_pcage_in(:)   ! pcnst-space
   integer,  intent(in) :: lspectoo_pcage_in(:)   ! pcnst-space
   real(r8), intent(in) :: mw_soa_host(:)
   real(r8), intent(in) :: mw_poa_host(:)
   real(r8), intent(in) :: rair               ! dry-air gas constant from host (J/K/kg)
   real(r8), intent(in) :: mwdry              ! dry-air molecular weight from host (kg/kmol)
   real(r8), intent(in) :: r_universal        ! universal gas constant from host (J/K/kmol)
   character(len=*), intent(out) :: errmsg
   integer,  intent(out) :: errflg

   ! local
   integer :: jsoa, l, l1, l2, n
   real(r8) :: tmp2

!-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   ! Store configuration
   ntot_amode_m = ntot_amode
   nsoa_m = nsoa
   npoa_m = npoa
   nspec_max_m = nspec_max

   ! Allocate and store mode configuration arrays
   allocate(nspec_amode_m(ntot_amode))
   nspec_amode_m(:) = nspec_amode(1:ntot_amode)

   allocate(alnsg_amode_m(ntot_amode))
   alnsg_amode_m(:) = alnsg_amode(1:ntot_amode)

   allocate(sigmag_amode_m(ntot_amode))
   sigmag_amode_m(:) = sigmag_amode(1:ntot_amode)

   allocate(specmw_amode_m(nspec_max, ntot_amode))
   specmw_amode_m(:,:) = specmw_amode(1:nspec_max, 1:ntot_amode)

   allocate(specdens_amode_m(nspec_max, ntot_amode))
   specdens_amode_m(:,:) = specdens_amode(1:nspec_max, 1:ntot_amode)

   ! Store species indices
   idx_h2so4_m = idx_h2so4
   idx_nh3_m = idx_nh3
   idx_msa_m = idx_msa

   allocate(idx_soag_m(nsoa))
   idx_soag_m(:) = idx_soag(1:nsoa)

   allocate(idx_so4_a_m(ntot_amode))
   idx_so4_a_m(:) = idx_so4_a(1:ntot_amode)

   allocate(idx_nh4_a_m(ntot_amode))
   idx_nh4_a_m(:) = idx_nh4_a(1:ntot_amode)

   allocate(idx_soa_a_m(ntot_amode, nsoa))
   idx_soa_a_m(:,:) = idx_soa_a(1:ntot_amode, 1:nsoa)

   allocate(idx_pom_a_m(ntot_amode, npoa))
   idx_pom_a_m(:,:) = idx_pom_a(1:ntot_amode, 1:npoa)

   allocate(idx_num_m(ntot_amode))
   idx_num_m(:) = idx_num(1:ntot_amode)

   allocate(idx_mass_m(nspec_max, ntot_amode))
   idx_mass_m(:,:) = idx_mass(1:nspec_max, 1:ntot_amode)

   ! Store pcarbon mode pointer
   modeptr_pcarbon_m = modeptr_pcarbon

   ! Store molecular weights
   allocate(mw_soa_host_m(nsoa))
   mw_soa_host_m(:) = mw_soa_host(1:nsoa)

   allocate(mw_poa_host_m(npoa))
   mw_poa_host_m(:) = mw_poa_host(1:npoa)

   ! Store host physical constants
   rair_m  = rair
   mwdry_m = mwdry
   rgas_m  = r_universal * 1.0e-3_r8  ! J/K/kmol -> J/K/mol

   ! Validate H2SO4 index (required species)
   if ((idx_h2so4 <= 0) .or. (idx_h2so4 > pcnst_in)) then
      write(errmsg, '(a,i7)') &
         'modal_aero_gasaerexch_init -- cannot find H2SO4 species, idx=', idx_h2so4
      errflg = 1
      return
   end if

   ! Compute species presence flags
   do_nh4g_m = .false.
   if ((idx_nh3 > 0) .and. (idx_nh3 <= pcnst_in)) do_nh4g_m = .true.

   do_msag_m = .false.
   if ((idx_msa > 0) .and. (idx_msa <= pcnst_in)) do_msag_m = .true.

   allocate(do_soag_m(nsoa))
   do_soag_any_m = .false.
   do_soag_m(:) = .false.
   do jsoa = 1, nsoa
      if ((idx_soag(jsoa) > 0) .and. (idx_soag(jsoa) <= pcnst_in)) then
         do_soag_any_m = .true.
         do_soag_m(jsoa) = .true.
      end if
   end do

  ! Compute ido arrays (species presence in modes)
   allocate(ido_so4a_m(ntot_amode))
   allocate(ido_nh4a_m(ntot_amode))
   allocate(ido_soaa_m(ntot_amode, nsoa))
   ido_so4a_m(:) = 0
   ido_nh4a_m(:) = 0
   ido_soaa_m(:,:) = 0

   ntot_soamode_m = 0
   do n = 1, ntot_amode
      l = idx_so4_a(n)
      if ((l > 0) .and. (l <= pcnst_in)) then
         ido_so4a_m(n) = 1
         if ( do_nh4g_m ) then
            l = idx_nh4_a(n)
            if ((l > 0) .and. (l <= pcnst_in)) then
               ido_nh4a_m(n) = 1
            end if
         end if
      end if

      do jsoa = 1, nsoa
         if ( do_soag_m(jsoa) ) then
            l = idx_soa_a(n,jsoa)
            if ((l > 0) .and. (l <= pcnst_in)) then
               ido_soaa_m(n,jsoa) = 1
               ntot_soamode_m = n
            end if
         end if
      end do ! jsoa
   end do ! n

  !
  !   define "from mode" and "to mode" for primary carbon aging
  !
  !   skip (turn off) aging if either is absent,
  !      or if accum mode so4 is absent
  !
   maxspec_pcage = nspec_max
   allocate(lspecfrm_pcage(maxspec_pcage))
   allocate(lspectoo_pcage(maxspec_pcage))
   allocate(soa_equivso4_factor(nsoa))
   allocate(fac_m2v_soa(nsoa))
   allocate(fac_m2v_pcarbon(nspec_max))

   lspecfrm_pcage(:) = 0
   lspectoo_pcage(:) = 0

   modefrm_pcage = -999888777
   modetoo_pcage = -999888777
   nspecfrm_pcage = 0

   if ((modeptr_pcarbon > 0) .and. (modeptr_accum > 0)) then
      l = idx_so4_a(modeptr_accum)
      if ((l > 0) .and. (l <= pcnst_in)) then
         modefrm_pcage = modeptr_pcarbon
         modetoo_pcage = modeptr_accum

         nspecfrm_pcage = nspecfrm_pcage_in
         lspecfrm_pcage(1:nspecfrm_pcage) = lspecfrm_pcage_in(1:nspecfrm_pcage)
         lspectoo_pcage(1:nspecfrm_pcage) = lspectoo_pcage_in(1:nspecfrm_pcage)
      end if
   end if

   if ( do_soag_any_m ) ntot_soamode_m = max( ntot_soamode_m, modefrm_pcage )

  ! Modify ido arrays for pcage mode
   if (modefrm_pcage > 0) then
      ido_so4a_m(modefrm_pcage) = 2
      if (ido_nh4a_m(modetoo_pcage) == 1) ido_nh4a_m(modefrm_pcage) = 2
      do jsoa = 1, nsoa
         if (ido_soaa_m(modetoo_pcage,jsoa) == 1) ido_soaa_m(modefrm_pcage,jsoa) = 2
      end do
   end if

  ! set for used in aging calcs:
  !    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
  !    soa_equivso4_factor(:)
   soa_equivso4_factor = 0.0_r8
   if (modefrm_pcage > 0) then
      n = modeptr_accum
      l2 = -1
      do l1 = 1, nspec_amode(n)
         if (idx_mass(l1,n) == idx_so4_a(n)) then
!               l2 = lspectype_amode(l1,n)
            l2 = l1
!               fac_m2v_so4 = specmw_amode(l2) / specdens_amode(l2)
            fac_m2v_so4 = specmw_amode(l1,n) / specdens_amode(l1,n)
!               tmp2 = spechygro(l2)
            tmp2 = spechygro(l1,n)

         end if
      end do
      if (l2 <= 0) then
         errmsg = 'modal_aero_gasaerexch_init error a002 finding accum. so4'
         errflg = 1
         return
      end if

      l2 = -1
      if (idx_nh4_a(n) > 0) then
         do l1 = 1, nspec_amode(n)
            if (idx_mass(l1,n) == idx_nh4_a(n)) then
!                  l2 = lspectype_amode(l1,n)
               l2 = l1
!                  fac_m2v_nh4 = specmw_amode(l2) / specdens_amode(l2)
               fac_m2v_nh4 = specmw_amode(l1,n) / specdens_amode(l1,n)

            end if
         end do
         if (l2 <= 0) then
            errmsg = 'modal_aero_gasaerexch_init error a002 finding accum. nh4'
            errflg = 1
            return
         end if
      else
         fac_m2v_nh4 = fac_m2v_so4
      end if

      do jsoa = 1, nsoa
         l2 = -1
         if (idx_soa_a(n,jsoa) <= 0) then
            write( errmsg, '(a,i4)') 'modal_aero_gasaerexch_init error a001 finding accum. jsoa =', jsoa
            errflg = 1
            return
         end if
         do l1 = 1, nspec_amode(n)
            if (idx_mass(l1,n) == idx_soa_a(n,jsoa)) then
!                  l2 = lspectype_amode(l1,n)
               l2 = l1
!                  fac_m2v_soa(jsoa) = specmw_amode(l2) / specdens_amode(l2)
               fac_m2v_soa(jsoa) = specmw_amode(l1,n) / specdens_amode(l1,n)
!                  soa_equivso4_factor(jsoa) = spechygro(l2)/tmp2
               soa_equivso4_factor(jsoa) = spechygro(l1,n)/tmp2
            end if
         end do
         if (l2 <= 0) then
            write( errmsg, '(a,i4)') 'modal_aero_gasaerexch_init error a002 finding accum. jsoa =', jsoa
            errflg = 1
            return
         end if
      end do

      fac_m2v_pcarbon(:) = 0.0_r8
      n = modeptr_pcarbon
      do l = 1, nspec_amode(n)
!            l2 = lspectype_amode(l,n)
!      fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!           [m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
!            fac_m2v_pcarbon(l) = specmw_amode(l2) / specdens_amode(l2)
         fac_m2v_pcarbon(l) = specmw_amode(l,n) / specdens_amode(l,n)
      end do
   end if

end subroutine modal_aero_gasaerexch_init

subroutine modal_aero_gasaerexch_run(                            &
                        ncol, pver, deltat, top_lev,             &
                        loffset,                                 &
                        t,        pmid,      pdel,     gravit,   &
                        troplev,                                 &
                        dgncur_a,           dgncur_awet,         &
                        use_sulfeq, sulfeq,                      &
                        num_q,                                   &
                        q,                                       &
                        dqdt, dotend, qsrflx_gaexch,             &
                        errmsg, errflg)
   integer,  intent(in)    :: ncol                 ! # of atmospheric columns
   integer,  intent(in)    :: pver                 ! # of vertical levels
   real(r8), intent(in)    :: deltat               ! time step [s]
   integer,  intent(in)    :: top_lev              ! top level for aerosol processes
   integer,  intent(in)    :: loffset              ! offset to convert pcnst-space to vmr-space [index]
   integer,  intent(in)    :: troplev(:)           ! (ncol) tropopause vertical index [index]
   real(r8), intent(in)    :: t(:,:)               ! (ncol,pver) temperature [K]
   real(r8), intent(in)    :: pmid(:,:)            ! (ncol,pver) pressure [Pa]
   real(r8), intent(in)    :: pdel(:,:)            ! (ncol,pver) pressure thickness [Pa]
   real(r8), intent(in)    :: gravit               ! gravitational acceleration [m s-2]
   real(r8), intent(in)    :: dgncur_a(:,:,:)      ! (ncol,pver,ntot_amode) dry diameter
   real(r8), intent(in)    :: dgncur_awet(:,:,:)   ! (ncol,pver,ntot_amode) wet diameter
   logical,  intent(in)    :: use_sulfeq           ! whether to use strat equilibrium
   real(r8), intent(in)    :: sulfeq(:,:,:)        ! (ncol,pver,ntot_amode) sulfeq values
   integer,  intent(in)    :: num_q                ! # of species in vmr array (= gas_pcnst)
   real(r8), intent(in)    :: q(:,:,:)             ! (ncol,pver,num_q) tracer VMR
                                                   ! *** MUST BE  #/kmol-air for number
                                                   ! *** MUST BE mol/mol-air for mass
   real(r8), intent(out)   :: dqdt(:,:,:)          ! (ncol,pver,num_q) tendencies
   logical,  intent(out)   :: dotend(:)            ! (num_q) which species have tendencies
   real(r8), intent(out)   :: qsrflx_gaexch(:,:)   ! (ncol,num_q) column-integrated gas-aerosol
                                                   ! exchange source/sink (kg/m2/s, pre adv_mass/mwdry
                                                   ! scaling) for the _sfgaex1 diagnostic.
                                                   ! Accumulated per-term here because the per-mode
                                                   ! and primary-carbon-aging contributions must be
                                                   ! summed separately to stay bfb with CAM.
   character(len=*), intent(out) :: errmsg
   integer,  intent(out)   :: errflg

   ! computes TMR (tracer mixing ratio) tendencies for gas condensation
   !    onto aerosol particles
   !
   ! this version does condensation of H2SO4, NH3, and MSA, both treated as
   ! completely non-volatile (gas --> aerosol, but no aerosol --> gas)
   !    gas H2SO4 goes to aerosol SO4
   !    gas MSA (if present) goes to aerosol SO4
   !       aerosol MSA is not distinguished from aerosol SO4
   !    gas NH3 (if present) goes to aerosol NH4
   !       if gas NH3 is not present, then ????

   ! local variables
   integer, parameter :: method_soa = 2
   !     method_soa=0 is no uptake
   !     method_soa=1 is irreversible uptake done like h2so4 uptake
   !     method_soa=2 is reversible uptake using subr modal_aero_soaexch

   integer :: i, iq, itmpa
   integer :: ido_so4a(ntot_amode_m), ido_nh4a(ntot_amode_m)
   integer ::  ido_soaa(ntot_amode_m,nsoa_m)
   integer :: j, jsoa
   integer :: k
   integer :: l, lsfrm, lstoo
   integer :: l_so4g, l_nh4g, l_msag
   integer :: l_soag(nsoa_m)
   integer :: n, niter, niter_max, ntot_soamode

   ! Local offset-adjusted index arrays (pcnst-space - loffset = vmr space)
   integer :: idx_so4_a_q(ntot_amode_m), idx_nh4_a_q(ntot_amode_m)
   integer :: idx_soa_a_q(ntot_amode_m,nsoa_m), idx_pom_a_q(ntot_amode_m,npoa_m)
   integer :: idx_num_q(ntot_amode_m), idx_mass_q(nspec_max_m,ntot_amode_m)
   integer :: lspecfrm_q(maxspec_pcage), lspectoo_q(maxspec_pcage)

   real (r8) :: avg_uprt_nh4, avg_uprt_so4, avg_uprt_soa(nsoa_m)
   real (r8) :: deltatxx
   real (r8) :: dqdt_nh4(ntot_amode_m), dqdt_so4(ntot_amode_m)
   real (r8) :: dqdt_soa(ntot_amode_m,nsoa_m)
   real (r8) :: dqdt_soag(nsoa_m)
   real (r8) :: fac_volsfc_pcarbon
   real (r8) :: fgain_nh4(ntot_amode_m), fgain_so4(ntot_amode_m)
   real (r8) :: fgain_soa(ntot_amode_m,nsoa_m)
   real(r8)  :: mw_poa_host(npoa_m)          ! molec wght of poa used in host code
   real(r8)  :: mw_soa_host(nsoa_m)          ! molec wght of poa used in host code
   real (r8) :: qmax_nh4, qnew_nh4, qnew_so4
   real (r8) :: qold_nh4(ntot_amode_m), qold_so4(ntot_amode_m)
   real (r8) :: qold_poa(ntot_amode_m,npoa_m)
   real (r8) :: qold_soa(ntot_amode_m,nsoa_m)
   real (r8) :: qold_soag(nsoa_m)
   real (r8) :: sum_dqdt_msa, sum_dqdt_so4
   real (r8) :: sum_dqdt_soa(nsoa_m)
   real (r8) :: sum_dqdt_nh4, sum_dqdt_nh4_b
   real (r8) :: sum_uprt_nh4, sum_uprt_so4
   real (r8) :: sum_uprt_soa(nsoa_m)
   real (r8) :: pdel_fac        ! pdel/gravit (kg/m2 per Pa), for column-integrated diagnostics
   real (r8) :: tmp1, tmp2, tmpa
   real (r8) :: tmp_kxt, tmp_pxt
   real (r8) :: tmp_so4a_bgn, tmp_so4a_end
   real (r8) :: tmp_so4g_avg, tmp_so4g_bgn, tmp_so4g_equ
   real (r8) :: uptkrate(ntot_amode_m,ncol,pver)
   real (r8) :: uptkratebb(ntot_amode_m)
   real (r8) :: uptkrate_soa(ntot_amode_m,nsoa_m)
                ! gas-to-aerosol mass transfer rates (1/s)
   real (r8) :: vol_core, vol_shell
   real (r8) :: xferfrac_pcage, xferfrac_max
   real (r8) :: xferrate

   logical  :: do_msag         ! true if msa gas is a species
   logical  :: do_nh4g         ! true if nh3 gas is a species
   logical  :: do_soag_any         ! true if soa gas is a species
   logical  :: do_soag(nsoa_m)       ! true if soa gas is a species


!----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

! set gas species indices from module-level storage, applying -loffset
! to convert pcnst-space to vmr (gas_pcnst) space
   l_so4g = idx_h2so4_m - loffset
   l_nh4g = idx_nh3_m - loffset
   l_msag = idx_msa_m - loffset
   do_nh4g = do_nh4g_m
   do_msag = do_msag_m
   do_soag_any = do_soag_any_m
   do_soag(:) = do_soag_m(:)
   do jsoa = 1, nsoa_m
      l_soag(jsoa) = idx_soag_m(jsoa) - loffset
   end do

! compute offset-adjusted per-mode index arrays
   idx_so4_a_q(:) = idx_so4_a_m(:) - loffset
   idx_nh4_a_q(:) = idx_nh4_a_m(:) - loffset
   idx_soa_a_q(:,:) = idx_soa_a_m(:,:) - loffset
   idx_pom_a_q(:,:) = idx_pom_a_m(:,:) - loffset
   idx_num_q(:) = idx_num_m(:) - loffset
   idx_mass_q(:,:) = idx_mass_m(:,:) - loffset
   do iq = 1, nspecfrm_pcage
      lspecfrm_q(iq) = lspecfrm_pcage(iq) - loffset
      lspectoo_q(iq) = lspectoo_pcage(iq)
      if (lspectoo_q(iq) > 0) lspectoo_q(iq) = lspectoo_q(iq) - loffset
   end do

! copy ido arrays from module-level storage
   ido_so4a(:) = ido_so4a_m(:)
   ido_nh4a(:) = ido_nh4a_m(:)
   ido_soaa(:,:) = ido_soaa_m(:,:)
   ntot_soamode = ntot_soamode_m

! set molecular weights from module-level storage
   mw_soa_host(:) = mw_soa_host_m(:)
   mw_poa_host(:) = mw_poa_host_m(:)

! set tendency flags
   dotend(:) = .false.

   dotend(l_so4g) = .true.
   if ( do_nh4g ) dotend(l_nh4g) = .true.
   if ( do_msag ) dotend(l_msag) = .true.
   do jsoa = 1, nsoa_m
      if ( do_soag(jsoa) ) dotend(l_soag(jsoa)) = .true.
   end do

   do n = 1, ntot_amode_m
      if (ido_so4a(n) == 1) then
         l = idx_so4_a_q(n)
         dotend(l) = .true.
         if ( do_nh4g ) then
            if (ido_nh4a(n) == 1) then
               l = idx_nh4_a_q(n)
               dotend(l) = .true.
            end if
         end if
      end if

      do jsoa = 1, nsoa_m
         if ( do_soag(jsoa) ) then
            if (ido_soaa(n,jsoa) == 1) then
               l = idx_soa_a_q(n,jsoa)
               dotend(l) = .true.
            end if
         end if
      end do ! jsoa
   end do ! n


   if (modefrm_pcage > 0) then
      do iq = 1, nspecfrm_pcage
         lsfrm = lspecfrm_q(iq)
         lstoo = lspectoo_q(iq)
         if ((lsfrm > 0) .and. (lsfrm <= num_q)) then
            dotend(lsfrm) = .true.
            if ((lstoo > 0) .and. (lstoo <= num_q)) then
               dotend(lstoo) = .true.
            end if
         end if
      end do


      n = modeptr_pcarbon_m
      fac_volsfc_pcarbon = exp( 2.5_r8*(alnsg_amode_m(n)**2) )
      xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps
   end if


! zero out tendencies
   dqdt(:,:,:) = 0.0_r8
   qsrflx_gaexch(:,:) = 0.0_r8

! compute gas-to-aerosol mass transfer rates
   call gas_aer_uptkrates( ncol,       pver,       top_lev,    &
                           loffset,                            &
                           q,          t,          pmid,       &
                           dgncur_awet,            uptkrate    )


! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)


   do k=top_lev,pver
      do i=1,ncol

!   fgain_so4(n) = fraction of total h2so4 uptake going to mode n
!   fgain_nh4(n) = fraction of total  nh3  uptake going to mode n
         sum_uprt_so4 = 0.0_r8
         sum_uprt_nh4 = 0.0_r8
         sum_uprt_soa = 0.0_r8
         do n = 1, ntot_amode_m
            uptkratebb(n) = uptkrate(n,i,k)
            if (ido_so4a(n) > 0) then
               fgain_so4(n) = uptkratebb(n)
               sum_uprt_so4 = sum_uprt_so4 + fgain_so4(n)
               if (ido_so4a(n) == 1) then
                  qold_so4(n) = q(i,k,idx_so4_a_q(n))
               else
                  qold_so4(n) = 0.0_r8
               end if
            else
               fgain_so4(n) = 0.0_r8
               qold_so4(n) = 0.0_r8
            end if

            if (ido_nh4a(n) > 0) then
               !   2.08 factor is for gas diffusivity (nh3/h2so4)
               !   differences in fuch-sutugin and accom coef ignored
               fgain_nh4(n) = uptkratebb(n)*2.08_r8
               sum_uprt_nh4 = sum_uprt_nh4 + fgain_nh4(n)
               if (ido_nh4a(n) == 1) then
                  qold_nh4(n) = q(i,k,idx_nh4_a_q(n))
               else
                  qold_nh4(n) = 0.0_r8
               end if
            else
               fgain_nh4(n) = 0.0_r8
               qold_nh4(n) = 0.0_r8
            end if

            do j = 1, npoa_m
               l = idx_pom_a_q(n,j)
               if (l > 0) then
                  qold_poa(n,j) = q(i,k,l)
               else
                  qold_poa(n,j) = 0.0_r8
               end if
            end do

            itmpa = 0
            do jsoa = 1, nsoa_m
               if (ido_soaa(n,jsoa) > 0) then
                  ! 0.81 factor is for gas diffusivity (soa/h2so4)
                  ! (differences in fuch-sutugin and accom coef ignored)
                  fgain_soa(n,jsoa) = uptkratebb(n)*0.81_r8
                  sum_uprt_soa(jsoa) = sum_uprt_soa(jsoa) + fgain_soa(n,jsoa)
                  if (ido_soaa(n,jsoa) == 1) then
                     l = idx_soa_a_q(n,jsoa)
                     qold_soa(n,jsoa) = q(i,k,l)
                     itmpa = itmpa + 1
                  else
                     qold_soa(n,jsoa) = 0.0_r8
                  end if
               else
                  fgain_soa(n,jsoa) = 0.0_r8
                  qold_soa(n,jsoa) = 0.0_r8
               end if
               uptkrate_soa(n,jsoa) = fgain_soa(n,jsoa)
            end do ! jsoa
            ! in previous code versions with nsoa=1,
            !    qold_poa was non-zero (i.e., loaded from q) only when ido_soaa(n)=1
            ! thus qold_poa=0 for the primary carbon mode which has ido_soaa=2
            ! this is probably not how it should be
            if (itmpa == 0) qold_poa(n,:) = 0.0_r8

         end do ! n

         if (sum_uprt_so4 > 0.0_r8) then
            do n = 1, ntot_amode_m
               fgain_so4(n) = fgain_so4(n) / sum_uprt_so4
            end do
         end if
!       at this point (sum_uprt_so4 <= 0.0) only when all the fgain_so4 are zero
         if (sum_uprt_nh4 > 0.0_r8) then
            do n = 1, ntot_amode_m
               fgain_nh4(n) = fgain_nh4(n) / sum_uprt_nh4
            end do
         end if

         do jsoa = 1, nsoa_m
            if (sum_uprt_soa(jsoa) > 0.0_r8) then
               do n = 1, ntot_amode_m
                  fgain_soa(n,jsoa) = fgain_soa(n,jsoa) / sum_uprt_soa(jsoa)
               end do
            end if
         end do

!   uptake amount (fraction of gas uptaken) over deltat
         avg_uprt_so4 = (1.0_r8 - exp(-deltatxx*sum_uprt_so4))/deltatxx
         avg_uprt_nh4 = (1.0_r8 - exp(-deltatxx*sum_uprt_nh4))/deltatxx

         do jsoa = 1, nsoa_m
            avg_uprt_soa(jsoa) = (1.0_r8 - exp(-deltatxx*sum_uprt_soa(jsoa)))/deltatxx
         end do

!   sum_dqdt_so4 = so4_a tendency from h2so4 gas uptake (mol/mol/s)
!   sum_dqdt_msa = msa_a tendency from msa   gas uptake (mol/mol/s)
!   sum_dqdt_nh4 = nh4_a tendency from nh3   gas uptake (mol/mol/s)
!   sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)
         sum_dqdt_so4 = q(i,k,l_so4g) * avg_uprt_so4
         if ( do_msag ) then
            sum_dqdt_msa = q(i,k,l_msag) * avg_uprt_so4
         else
            sum_dqdt_msa = 0.0_r8
         end if
         if ( do_nh4g ) then
            sum_dqdt_nh4 = q(i,k,l_nh4g) * avg_uprt_nh4
         else
            sum_dqdt_nh4 = 0.0_r8
         end if

         do jsoa = 1, nsoa_m
            if ( do_soag(jsoa) ) then
               sum_dqdt_soa(jsoa) = q(i,k,l_soag(jsoa)) * avg_uprt_soa(jsoa)
            else
               sum_dqdt_soa(jsoa) = 0.0_r8
            end if
         end do

         if ( use_sulfeq .and. (k <= troplev(i)) ) then
            !   compute TMR tendencies for so4 interstial aerosol due to reversible gas uptake
            !   only above the tropopause

            tmp_kxt = deltatxx*sum_uprt_so4  ! sum over modes of uptake_rate*deltat
            tmp_pxt = 0.0_r8
            do n = 1, ntot_amode_m
               if (ido_so4a(n) <= 0) cycle
               tmp_pxt = tmp_pxt + uptkratebb(n)*sulfeq(i,k,n)
            end do
            tmp_pxt = max( 0.0_r8, tmp_pxt*deltatxx )  ! sum over modes of uptake_rate*sulfeq*deltat
            tmp_so4g_bgn = q(i,k,l_so4g)
            ! calc avg h2so4(g) over deltat
            if (tmp_kxt >= 1.0e-5_r8) then
               ! exponential decay towards equilibrium value solution
               tmp_so4g_equ = tmp_pxt/tmp_kxt
               tmp_so4g_avg = tmp_so4g_equ + (tmp_so4g_bgn-tmp_so4g_equ)*(1.0_r8-exp(-tmp_kxt))/tmp_kxt
            else
               ! first order approx for tmp_kxt small
               tmp_so4g_avg = tmp_so4g_bgn*(1.0_r8-0.5_r8*tmp_kxt) + 0.5_r8*tmp_pxt
            end if
            sum_dqdt_so4 = 0.0_r8
            do n = 1, ntot_amode_m
               if (ido_so4a(n) <= 0) cycle
               ! calc change to so4(a) in mode n
               if (ido_so4a(n) == 1) then
                  l = idx_so4_a_q(n)
                  tmp_so4a_bgn = q(i,k,l)
               else
                  tmp_so4a_bgn = 0.0_r8
               end if
               tmp_so4a_end = tmp_so4a_bgn + deltatxx*uptkratebb(n)*(tmp_so4g_avg-sulfeq(i,k,n))
               tmp_so4a_end = max( 0.0_r8, tmp_so4a_end )
               dqdt_so4(n) = (tmp_so4a_end - tmp_so4a_bgn)/deltatxx
               sum_dqdt_so4 = sum_dqdt_so4 + dqdt_so4(n)
            end do
            ! do not allow msa condensation in stratosphere
            ! ( Note that the code for msa has never been used.
            !   The plan was to simulate msa(g), treat it as non-volatile (like h2so4(g)),
            !   and treat condensed msa as sulfate, so just one additional tracer. )
            if ( do_msag ) sum_dqdt_msa = 0.0_r8

         else
            !   compute TMR tendencies for so4 interstial aerosol due to simple gas uptake
            do n = 1, ntot_amode_m
               dqdt_so4(n) = fgain_so4(n)*(sum_dqdt_so4 + sum_dqdt_msa)
            end do
         end if

         !   compute TMR tendencies for nh4 interstial aerosol due to simple gas uptake
         !   but force nh4/so4 molar ratio <= 2
         sum_dqdt_nh4_b = 0.0_r8
         dqdt_nh4(:) = 0._r8
         if ( do_nh4g ) then
            do n = 1, ntot_amode_m
               dqdt_nh4(n) = fgain_nh4(n)*sum_dqdt_nh4
               qnew_nh4 = qold_nh4(n) + dqdt_nh4(n)*deltat
               qnew_so4 = qold_so4(n) + dqdt_so4(n)*deltat
               qmax_nh4 = 2.0_r8*qnew_so4
               if (qnew_nh4 > qmax_nh4) then
                  dqdt_nh4(n) = (qmax_nh4 - qold_nh4(n))/deltatxx
               end if
               sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4(n)
            end do
         end if

         if (( do_soag_any ) .and. (method_soa > 1)) then
!   compute TMR tendencies for soag and soa interstial aerosol
!   using soa parameterization
            niter_max = 1000
            dqdt_soa(:,:) = 0.0_r8
            dqdt_soag(:) = 0.0_r8
            do jsoa = 1, nsoa_m
               qold_soag(jsoa) = q(i,k,l_soag(jsoa))
            end do

            call modal_aero_soaexch( deltat, t(i,k), pmid(i,k), &
                 niter, niter_max, ntot_amode_m, ntot_soamode, npoa_m, nsoa_m, &
                 mw_poa_host, mw_soa_host, &
                 qold_soag, qold_soa, qold_poa, uptkrate_soa, &
                 dqdt_soag, dqdt_soa )
            sum_dqdt_soa(:) = -dqdt_soag(:)

         else if ( do_soag_any ) then
!   compute TMR tendencies for soa interstial aerosol
!   due to simple gas uptake

            do jsoa = 1, nsoa_m
               do n = 1, ntot_amode_m
                  dqdt_soa(n,jsoa) = fgain_soa(n,jsoa)*sum_dqdt_soa(jsoa)
               end do
            end do
         else ! method_soa is neither 1 nor 2, no uptake
            dqdt_soa(:,:) = 0.0_r8
         end if

         pdel_fac = pdel(i,k)/gravit
         do n = 1, ntot_amode_m
            if (ido_so4a(n) == 1) then
               l = idx_so4_a_q(n)
               dqdt(i,k,l) = dqdt_so4(n)
               qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt_so4(n)*pdel_fac
            end if

            if ( do_nh4g ) then
               if (ido_nh4a(n) == 1) then
                  l = idx_nh4_a_q(n)
                  dqdt(i,k,l) = dqdt_nh4(n)
                  qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt_nh4(n)*pdel_fac
               end if
            end if

            do jsoa = 1, nsoa_m
               if ( do_soag(jsoa) ) then
                  if (ido_soaa(n,jsoa) == 1) then
                     l = idx_soa_a_q(n,jsoa)
                     dqdt(i,k,l) = dqdt_soa(n,jsoa) !calculated by  modal_aero_soaexch for method_soa=2
                     qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt_soa(n,jsoa)*pdel_fac
                  end if
               end if
            end do
         end do ! n

!   compute TMR tendencies for h2so4, nh3, and msa gas
!   due to simple gas uptake
         l = l_so4g
         dqdt(i,k,l) = -sum_dqdt_so4
         qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt(i,k,l)*pdel_fac

         if ( do_msag ) then
            l = l_msag
            dqdt(i,k,l) = -sum_dqdt_msa
            qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt(i,k,l)*pdel_fac
         end if

         if ( do_nh4g ) then
            l = l_nh4g
            dqdt(i,k,l) = -sum_dqdt_nh4_b
            qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt(i,k,l)*pdel_fac
         end if

         do jsoa = 1, nsoa_m
            if ( do_soag(jsoa) ) then
               l = l_soag(jsoa)
               dqdt(i,k,l) = -sum_dqdt_soa(jsoa)
! dqdt for gas is negative of the sum of dqdt for aerosol soa species in each mode: Manish
               qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt(i,k,l)*pdel_fac
            end if
         end do

!   compute TMR tendencies associated with primary carbon aging
         if (modefrm_pcage > 0) then
            n = modeptr_pcarbon_m
            tmpa = 0.0_r8
            do jsoa = 1, nsoa_m
               tmpa = tmpa + dqdt_soa(n,jsoa)*fac_m2v_soa(jsoa)*soa_equivso4_factor(jsoa)
            end do
            vol_shell = deltat *   &
                 ( dqdt_so4(n)*fac_m2v_so4 + dqdt_nh4(n)*fac_m2v_nh4 + tmpa )
            vol_core = 0.0_r8
            do l = 1, nspec_amode_m(n)
               vol_core = vol_core + &
                    q(i,k,idx_mass_q(l,n))*fac_m2v_pcarbon(l)
            end do
!   ratio1 = vol_shell/vol_core =
!      actual hygroscopic-shell-volume/carbon-core-volume after gas uptake
!   ratio2 = 6.0_r8*dr_so4_monolayers_pcage/(dgncur_a*fac_volsfc_pcarbon)
!      = (shell-volume corresponding to n_so4_monolayers_pcage)/core-volume
!      The 6.0/(dgncur_a*fac_volsfc_pcarbon) = (mode-surface-area/mode-volume)
!   Note that vol_shell includes both so4+nh4 AND soa as "equivalent so4",
!      The soa_equivso4_factor accounts for the lower hygroscopicity of soa.
!
!   Define xferfrac_pcage = min( 1.0, ratio1/ratio2)
!   But ratio1/ratio2 == tmp1/tmp2, and coding below avoids possible overflow
!
            tmp1 = vol_shell*dgncur_a(i,k,n)*fac_volsfc_pcarbon
            tmp2 = max( 6.0_r8*dr_so4_monolayers_pcage*vol_core, 0.0_r8 )
            if (tmp1 >= tmp2) then
               xferfrac_pcage = xferfrac_max
            else
               xferfrac_pcage = min( tmp1/tmp2, xferfrac_max )
            end if

            if (xferfrac_pcage > 0.0_r8) then
               do iq = 1, nspecfrm_pcage
                  lsfrm = lspecfrm_q(iq)
                  lstoo = lspectoo_q(iq)
                  xferrate = (xferfrac_pcage/deltat)*q(i,k,lsfrm)
                  dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferrate
                  qsrflx_gaexch(i,lsfrm) = qsrflx_gaexch(i,lsfrm) - xferrate*pdel_fac
                  if ((lstoo > 0) .and. (lstoo <= num_q)) then
                     dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferrate
                     qsrflx_gaexch(i,lstoo) = qsrflx_gaexch(i,lstoo) + xferrate*pdel_fac
                  end if
               end do

               if (ido_so4a(modetoo_pcage) > 0) then
                  l = idx_so4_a_q(modetoo_pcage)
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_so4(modefrm_pcage)
                  qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt_so4(modefrm_pcage)*pdel_fac
               end if

               if (ido_nh4a(modetoo_pcage) > 0) then
                  l = idx_nh4_a_q(modetoo_pcage)
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_nh4(modefrm_pcage)
                  qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt_nh4(modefrm_pcage)*pdel_fac
               end if

               do jsoa = 1, nsoa_m
                  if (ido_soaa(modetoo_pcage,jsoa) > 0) then
                     l = idx_soa_a_q(modetoo_pcage,jsoa)
                     dqdt(i,k,l) = dqdt(i,k,l) + dqdt_soa(modefrm_pcage,jsoa)
                     qsrflx_gaexch(i,l) = qsrflx_gaexch(i,l) + dqdt_soa(modefrm_pcage,jsoa)*pdel_fac
                  end if
               end do

            end if

         end if

      end do   ! "i = 1, ncol"
   end do     ! "k = top_lev, pver"

end subroutine modal_aero_gasaerexch_run


subroutine gas_aer_uptkrates( ncol,       pver,       top_lev,    &
                              loffset,                            &
                              q,          t,          pmid,       &
                              dgncur_awet,            uptkrate    )
   !
   !                         /
   !   computes   uptkrate = | dx  dN/dx  gas_conden_rate(Dp(x))
   !                         /
   !   using Gauss-Hermite quadrature of order nghq=2
   !
   !       Dp = particle diameter (cm)
   !       x = ln(Dp)
   !       dN/dx = log-normal particle number density distribution
   !       gas_conden_rate(Dp) = 2 * pi * gasdiffus * Dp * F(Kn,ac)
   !           F(Kn,ac) = Fuchs-Sutugin correction factor
   !           Kn = Knudsen number
   !           ac = accomodation coefficient
   !

   integer,  intent(in) :: ncol                 ! number of atmospheric column
   integer,  intent(in) :: pver                 ! number of vertical levels
   integer,  intent(in) :: top_lev              ! top level for aerosol processes
   integer,  intent(in) :: loffset              ! offset to convert pcnst-space to vmr space
   real(r8), intent(in) :: q(:,:,:)             ! (ncol,pver,num_q) Tracer array (mol,#/mol-air)
   real(r8), intent(in) :: t(:,:)               ! (ncol,pver) Temperature in Kelvin
   real(r8), intent(in) :: pmid(:,:)            ! (ncol,pver) Air pressure in Pa
   real(r8), intent(in) :: dgncur_awet(:,:,:)   ! (ncol,pver,ntot_amode_m)

   real(r8), intent(out) :: uptkrate(:,:,:)     ! (ntot_amode_m,ncol,pver)
                            ! gas-to-aerosol mass transfer rates (1/s)


! local
   integer, parameter :: nghq = 2
   integer :: i, iq, k, n

   ! Can use sqrt here once Lahey is gone.
   real(r8), parameter :: tworootpi = 3.5449077_r8
   real(r8), parameter :: root2 = 1.4142135_r8
   real(r8), parameter :: beta = 2.0_r8

   real(r8) :: aircon
   real(r8) :: const
   real(r8) :: dp
   real(r8) :: gasdiffus, gasspeed
   real(r8) :: freepathx2, fuchs_sutugin
   real(r8) :: knudsen
   real(r8) :: lndp, lndpgn, lnsg
   real(r8) :: num_a
   real(r8) :: rhoair
   real(r8) :: sumghq
   real(r8), save :: xghq(nghq), wghq(nghq) ! quadrature abscissae and weights

   data xghq / 0.70710678_r8, -0.70710678_r8 /
   data wghq / 0.88622693_r8,  0.88622693_r8 /


! outermost loop over all modes
   do n = 1, ntot_amode_m

! 22-aug-2007 rc easter - get number from q array rather
!    than computing a "bounded" number conc.
!! compute dry volume = sum_over_components{ component_mass / density }
!!    (m3-AP/mol-air)
!! compute it for all i,k to improve accessing q array
!      dryvol_a(1:ncol,:) = 0.0_r8
!      do l1 = 1, nspec_amode(n)
!         l2 = lspectype_amode(l1,n)
!! dum_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!! [m3-AP/kmol-AP]= [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
!         dum_m2v = specmw_amode(l2) / specdens_amode(l2)
!         la = lmassptr_amode(l1,n)
!         dryvol_a(1:ncol,:) = dryvol_a(1:ncol,:)    &
!                            + max(0.0_r8,q(1:ncol,:,la))*dum_m2v
!      end do

! loops k and i
      do k=top_lev,pver
      do i=1,ncol

         rhoair = pmid(i,k)/(rair_m*t(i,k))   ! (kg-air/m3)
!        aircon = 1.0e3*rhoair/mwdry        ! (mol-air/m3)

!!   "bounded" number conc. (#/m3)
!        num_a = dryvol_a(i,k)*v2ncur_a(i,k,n)*aircon

!   number conc. (#/m3) -- note q(i,k,numptr) is (#/kmol-air)
!   so need aircon in (kmol-air/m3)
         aircon = rhoair/mwdry_m            ! (kmol-air/m3)
         num_a = q(i,k,idx_num_m(n)-loffset)*aircon

!   gasdiffus = h2so4 gas diffusivity from mosaic code (m^2/s)
!               (pmid must be Pa)
         gasdiffus = 0.557e-4_r8 * (t(i,k)**1.75_r8) / pmid(i,k)
!   gasspeed = h2so4 gas mean molecular speed from mosaic code (m/s)
         gasspeed  = 1.470e1_r8 * sqrt(t(i,k))
!   freepathx2 = 2 * (h2so4 mean free path)  (m)
         freepathx2 = 6.0_r8*gasdiffus/gasspeed

         lnsg   = log( sigmag_amode_m(n) )
         lndpgn = log( dgncur_awet(i,k,n) )   ! (m)
         const  = tworootpi * num_a * exp(beta*lndpgn + 0.5_r8*(beta*lnsg)**2)

!   sum over gauss-hermite quadrature points
         sumghq = 0.0_r8
         do iq = 1, nghq
            lndp = lndpgn + beta*lnsg**2 + root2*lnsg*xghq(iq)
            dp = exp(lndp)

!   knudsen number
            knudsen = freepathx2/dp
!  Changed by Manish Shrivastava on 7/17/2013 to use accom=1; because we do not know better
!   following assumes accomodation coefficient = ac = 1. instead 0.65 ! answer change needs to be tested
!   (Adams & Seinfeld, 2002, JGR, and references therein)
!           fuchs_sutugin = (0.75*ac*(1. + knudsen)) /
!                           (knudsen*(1.0 + knudsen + 0.283*ac) + 0.75*ac)
            fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) /   &
                            (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)
            sumghq = sumghq + wghq(iq)*dp*fuchs_sutugin/(dp**beta)
         end do
         uptkrate(n,i,k) = const * gasdiffus * sumghq

      end do   ! "do i = 1, ncol"
      end do   ! "do k = 1, pver"

   end do   ! "do n = 1, ntot_soamode"

end subroutine gas_aer_uptkrates

!----------------------------------------------------------------------

subroutine modal_aero_soaexch( dtfull, temp, pres, &
          niter, niter_max, ntot_amode, ntot_soamode, ntot_poaspec, ntot_soaspec, &
          mw_poa_host, mw_soa_host, &
          g_soa_in, a_soa_in, a_poa_in, xferrate_in, &
          g_soa_tend, a_soa_tend )
!         g_soa_tend, a_soa_tend, g0_soa, idiagss )

!-----------------------------------------------------------------------
!
! Purpose:
!
! calculates condensation/evaporation of "soa gas"
! to/from multiple aerosol modes in 1 grid cell
!
! key assumptions
! (1) ambient equilibrium vapor pressure of soa gas
!     is given by p0_soa_298 and delh_vap_soa
! (2) equilibrium vapor pressure of soa gas at aerosol
!     particle surface is given by raoults law in the form
!     g_star = g0_soa*[a_soa/(a_soa + a_opoa)]
! (3) (oxidized poa)/(total poa) is equal to frac_opoa (constant)
!
!
! Author: R. Easter and R. Zaveri
! Additions to run with multiple BC, SOA and POM's: Shrivastava et al., 2015
!-----------------------------------------------------------------------

      implicit none

      real(r8), intent(in)  :: dtfull                     ! full integration time step (s)
      real(r8), intent(in)  :: temp                       ! air temperature (K)
      real(r8), intent(in)  :: pres                       ! air pressure (Pa)
      integer,  intent(out) :: niter                      ! number of iterations performed
      integer,  intent(in)  :: niter_max                  ! max allowed number of iterations
      integer,  intent(in)  :: ntot_amode                 ! number of modes
      integer,  intent(in)  :: ntot_soamode               ! number of modes having soa
      integer,  intent(in)  :: ntot_poaspec               ! number of poa species
      integer,  intent(in)  :: ntot_soaspec               ! number of soa species
      real(r8), intent(in)  :: mw_poa_host(ntot_poaspec)  ! molec wght of poa used in host code
      real(r8), intent(in)  :: mw_soa_host(ntot_soaspec)  ! molec wght of poa used in host code
      real(r8), intent(in)  :: g_soa_in(ntot_soaspec)               ! initial soa gas mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: a_soa_in(ntot_amode,ntot_soaspec)    ! initial soa aerosol mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: a_poa_in(ntot_amode,ntot_poaspec)    ! initial poa aerosol mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: xferrate_in(ntot_amode,ntot_soaspec) ! gas-aerosol mass transfer rate (1/s)
      real(r8), intent(out) :: g_soa_tend(ntot_soaspec)             ! soa gas mixrat tendency (mol/mol/s at host mw)
      real(r8), intent(out) :: a_soa_tend(ntot_amode,ntot_soaspec)  ! soa aerosol mixrat tendency (mol/mol/s at host mw)
!     integer,  intent(in)  :: idiagss

      integer :: ll
      integer :: m,k

      logical :: skip_soamode(ntot_amode)   ! true if this mode does not have soa

      real(r8) :: rgas
      real(r8), parameter :: a_min1 = 1.0e-20_r8
      real(r8), parameter :: g_min1 = 1.0e-20_r8
      real(r8), parameter :: alpha = 0.05_r8     ! parameter used in calc of time step
      real(r8), parameter :: dtsub_fixed = -1.0_r8  ! fixed sub-step for time integration (s)

      real(r8) :: a_ooa_sum_tmp(ntot_soamode)          ! total ooa (=soa+opoa) in a mode
      real(r8) :: a_opoa(ntot_soamode)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa(ntot_soamode,ntot_soaspec)     ! soa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa_tmp(ntot_soamode,ntot_soaspec) ! temporary soa aerosol mixrat (mol/mol)
      real(r8) :: beta(ntot_soamode,ntot_soaspec)      ! dtcur*xferrate
      real(r8) :: delh_vap_soa(ntot_soaspec)           ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
      real(r8) :: del_g_soa_tmp(ntot_soaspec)
      real(r8) :: dtcur                                ! current time step (s)
      real(r8) :: dtmax                                ! = (dtfull-tcur)
      real(r8) :: g0_soa(ntot_soaspec)                 ! ambient soa gas equilib mixrat (mol/mol at actual mw)
      real(r8) :: g_soa(ntot_soaspec)                  ! soa gas mixrat (mol/mol at actual mw)
      real(r8) :: g_star(ntot_soamode,ntot_soaspec)    ! soa gas mixrat that is in equilib
                                                       ! with each aerosol mode (mol/mol)
      real(r8) :: mw_poa(ntot_poaspec)                 ! actual molec wght of poa
      real(r8) :: mw_soa(ntot_soaspec)                 ! actual molec wght of soa
      real(r8) :: opoa_frac(ntot_poaspec)              ! fraction of poa that is opoa
      real(r8) :: phi(ntot_soamode,ntot_soaspec)       ! "relative driving force"
      real(r8) :: p0_soa(ntot_soaspec)                 ! soa gas equilib vapor presssure (atm)
      real(r8) :: p0_soa_298(ntot_soaspec)             ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k
      real(r8) :: sat(ntot_soamode,ntot_soaspec)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                       !    used by the numerical integration scheme -- it is not a saturation rato!
      real(r8) :: tcur                                 ! current integration time (from 0 s)
      real(r8) :: tmpa, tmpb, tmpf
      real(r8) :: tot_soa(ntot_soaspec)                ! g_soa + sum( a_soa(:) )
      real(r8) :: xferrate(ntot_amode,ntot_soaspec)    ! gas-aerosol mass transfer rate (1/s)

      rgas = rgas_m

      ! Changed by Manish Shrivastava
      opoa_frac(:) = 0.0_r8 !POA does not form solution with SOA for all runs; set opoa_frac=0.0_r8  by Manish Shrivastava
      mw_poa(:) = 250.0_r8
      mw_soa(:) = 250.0_r8

      ! New SOA properties added by Manish Shrivastava on 09/27/2012
      if (ntot_soaspec ==1) then
         p0_soa_298(:) = 9.7831E-11_r8
         delh_vap_soa(:) = 131.0e3_r8
         opoa_frac(:) = 0.0_r8
      elseif (ntot_soaspec ==2) then
         ! same for anthropogenic and biomass burning species
         p0_soa_298 (1) = 1.0e-10_r8
         p0_soa_298 (2) = 1.0e-10_r8
         delh_vap_soa(:) = 156.0e3_r8
      elseif(ntot_soaspec ==5) then
         ! 5 volatility bins for each of the a combined SOA classes ( including biomass burning, fossil fuel, biogenic)
         p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
         p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
         p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
         p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
         p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3

         delh_vap_soa(1) = 153.0e3_r8
         delh_vap_soa(2) = 142.0e3_r8
         delh_vap_soa(3) = 131.0e3_r8
         delh_vap_soa(4) = 120.0e3_r8
         delh_vap_soa(5) = 109.0e3_r8
      elseif(ntot_soaspec ==15) then
         !
         ! 5 volatility bins for each of the 3 SOA classes ( biomass burning, fossil fuel, biogenic)
         ! SOA species 1-5 are for anthropogenic while 6-10 are for biomass burning SOA
         ! SOA species 11-15 are for biogenic SOA, based on Cappa et al., Reference needs to be updated
         ! For MW=250.0
         p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
         p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
         p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
         p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
         p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3
         p0_soa_298 (6) = 9.7831E-13_r8 !soabb0 C*=0.01ug/m3
         p0_soa_298 (7) = 9.7831E-12_r8 !soabb1 C*=0.10ug/m3
         p0_soa_298 (8) = 9.7831E-11_r8 !soabb2 C*=1.0ug/m3
         p0_soa_298 (9) = 9.7831E-10_r8 !soabb3 C*=10.0ug/m3
         p0_soa_298 (10) = 9.7831E-9_r8  !soabb4 C*=100.0ug/m3
         p0_soa_298 (11) = 9.7831E-13_r8 !soabg0 C*=0.01ug/m3
         p0_soa_298 (12) = 9.7831E-12_r8 !soabg1 C*=0.1ug/m3
         p0_soa_298 (13) = 9.7831E-11_r8 !soabg2 C*=1.0ug/m3
         p0_soa_298 (14) = 9.7831E-10_r8 !soabg3 C*=10.0ug/m3
         p0_soa_298 (15) = 9.7831E-9_r8  !soabg4 C*=100.0ug/m3

         !
         ! have to be adjusted to 15 species, following the numbers by Epstein et al., 2012
         !
         delh_vap_soa(1) = 153.0e3_r8
         delh_vap_soa(2) = 142.0e3_r8
         delh_vap_soa(3) = 131.0e3_r8
         delh_vap_soa(4) = 120.0e3_r8
         delh_vap_soa(5) = 109.0e3_r8
         delh_vap_soa(6) = 153.0e3_r8
         delh_vap_soa(7) = 142.0e3_r8
         delh_vap_soa(8) = 131.0e3_r8
         delh_vap_soa(9) = 120.0e3_r8
         delh_vap_soa(10) = 109.0e3_r8
         delh_vap_soa(11) = 153.0e3_r8
         delh_vap_soa(12) = 142.0e3_r8
         delh_vap_soa(13) = 131.0e3_r8
         delh_vap_soa(14) = 120.0e3_r8
         delh_vap_soa(15) = 109.0e3_r8
      endif

      !BSINGH - Initialized g_soa_tend and a_soa_tend to circumvent the undefined behavior (04/16/12)
      g_soa_tend(:)   = 0.0_r8
      a_soa_tend(:,:) = 0.0_r8

      ! determine which modes have non-zero transfer rates
      !    and are involved in the soa gas-aerosol transfer
      ! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
      do m = 1, ntot_soamode
         skip_soamode(m) = .true.
         do ll = 1, ntot_soaspec
            xferrate(m,ll) = xferrate_in(m,ll)
            skip_soamode(m) = .false.
         end do
      end do

      ! convert incoming mixing ratios from mol/mol at the "host-code" molec. weight (12.0 in cam5)
      !    to mol/mol at the "actual" molec. weight (currently assumed to be 250.0)
      ! also
      !    force things to be non-negative
      !    calc tot_soa(ll)
      !    calc a_opoa (always slightly >0)
      do ll = 1, ntot_soaspec
         tmpf = mw_soa_host(ll)/mw_soa(ll)
         g_soa(ll) = max( g_soa_in(ll), 0.0_r8 ) * tmpf
         tot_soa(ll) = g_soa(ll)
         do m = 1, ntot_soamode
            if ( skip_soamode(m) ) cycle
            a_soa(m,ll) = max( a_soa_in(m,ll), 0.0_r8 ) * tmpf
            tot_soa(ll) = tot_soa(ll) + a_soa(m,ll)
         end do
      end do

      do m = 1, ntot_soamode
         if ( skip_soamode(m) ) cycle
         a_opoa(m) = 0.0_r8
         do ll = 1, ntot_poaspec
            a_opoa(m) = opoa_frac(ll)*a_poa_in(m,ll)
         end do
      end do

      ! calc ambient equilibrium soa gas
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
              exp( -(delh_vap_soa(ll)/rgas)*((1.0_r8/temp)-(1.0_r8/298.0_r8)) )
         g0_soa(ll) = 1.01325e5_r8*p0_soa(ll)/pres
      end do

      niter = 0
      tcur = 0.0_r8
      dtcur = 0.0_r8
      phi(:,:) = 0.0_r8
      g_star(:,:) = 0.0_r8

!     if (idiagss > 0) then
!        write(luna,'(a,1p,10e11.3)') 'p0, g0_soa', p0_soa, g0_soa
!        write(luna,'(3a)') &
!           'niter, tcur,   dtcur,    phi(:),                       ', &
!           'g_star(:),                    ', &
!           'a_soa(:),                     g_soa'
!        write(luna,'(3a)') &
!           '                         sat(:),                       ', &
!           'sat(:)*a_soa(:)               ', &
!           'a_opoa(:)'
!        write(luna,'(i3,1p,20e10.2)') niter, tcur, dtcur, &
!           phi(:), g_star(:), a_soa(:), g_soa
!     end if


! integration loop -- does multiple substeps to reach dtfull
time_loop: &
      do while (tcur < dtfull-1.0e-3_r8 )

      niter = niter + 1
      if (niter > niter_max) exit

      tmpa = 0.0_r8  ! time integration parameter for all soa species
      do m = 1, ntot_soamode
         if ( skip_soamode(m) ) cycle
         a_ooa_sum_tmp(m) = a_opoa(m) + sum( a_soa(m,1:ntot_soaspec) )
      end do
      do ll = 1, ntot_soaspec
         tmpb = 0.0_r8  ! time integration parameter for a single soa species
         do m = 1, ntot_soamode
            if ( skip_soamode(m) ) cycle
            sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
            g_star(m,ll) = sat(m,ll)*a_soa(m,ll)
            phi(m,ll) = (g_soa(ll) - g_star(m,ll))/max( g_soa(ll), g_star(m,ll), g_min1 )
            tmpb = tmpb + xferrate(m,ll)*abs(phi(m,ll))
         end do
         tmpa = max( tmpa, tmpb )
      end do

      if (dtsub_fixed > 0.0_r8) then
         dtcur = dtsub_fixed
         tcur = tcur + dtcur
      else
         dtmax = dtfull-tcur
         if (dtmax*tmpa <= alpha) then
! here alpha/tmpa >= dtmax, so this is final substep
            dtcur = dtmax
            tcur = dtfull
         else
            dtcur = alpha/tmpa
            tcur = tcur + dtcur
         end if
      end if

! step 1 - for modes where soa is condensing, estimate "new" a_soa(m,ll)
!    using an explicit calculation with "old" g_soa
!    and g_star(m,ll) calculated using "old" a_soa(m,ll)
! do this to get better estimate of "new" a_soa(m,ll) and sat(m,ll)
      do m = 1, ntot_soamode
         if ( skip_soamode(m) ) cycle
         do ll = 1, ntot_soaspec
            ! first ll loop calcs a_soa_tmp(m,ll) & a_ooa_sum_tmp
            a_soa_tmp(m,ll) = a_soa(m,ll)
            beta(m,ll) = dtcur*xferrate(m,ll)
            del_g_soa_tmp(ll) = g_soa(ll) - g_star(m,ll)
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               a_soa_tmp(m,ll) = a_soa(m,ll) + beta(m,ll)*del_g_soa_tmp(ll)
            end if
         end do
         a_ooa_sum_tmp(m) = a_opoa(m) + sum( a_soa_tmp(m,1:ntot_soaspec) )
         do ll = 1, ntot_soaspec
            ! second ll loop calcs sat & g_star
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
               g_star(m,ll) = sat(m,ll)*a_soa_tmp(m,ll)   ! this just needed for diagnostics
            end if
         end do
      end do

! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(m,ll) calculated semi-implicitly
      do ll = 1, ntot_soaspec
         tmpa = 0.0_r8
         tmpb = 0.0_r8
         do m = 1, ntot_soamode
            if ( skip_soamode(m) ) cycle
            tmpa = tmpa + a_soa(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
            tmpb = tmpb + beta(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
         end do

         g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_r8 + tmpb)
         g_soa(ll) = max( 0.0_r8, g_soa(ll) )
         do m = 1, ntot_soamode
            if ( skip_soamode(m) ) cycle
            a_soa(m,ll) = (a_soa(m,ll) + beta(m,ll)*g_soa(ll))/   &
                       (1.0_r8 + beta(m,ll)*sat(m,ll))
         end do
      end do

!     if (idiagss > 0) then
!        write(luna,'(i3,1p,20e10.2)') niter, tcur, dtcur, &
!           phi(:), g_star(:), a_soa(:), g_soa
!        write(luna,'(23x,1p,20e10.2)') &
!           sat(:), sat(:)*a_soa(:), a_opoa(:)
!     end if

!     if (niter > 9992000) then
!        write(luna,'(a)') '*** to many iterations'
!        exit
!     end if

      end do time_loop


! calculate outgoing tendencies (at the host-code molec. weight)
! (a_soa & g_soa are at actual mw, but a_soa_in & g_soa_in are at host-code mw)
      do ll = 1, ntot_soaspec
         tmpf = mw_soa(ll)/mw_soa_host(ll)
         g_soa_tend(ll) = (g_soa(ll)*tmpf - g_soa_in(ll))/dtfull
         do m = 1, ntot_soamode
            if ( skip_soamode(m) ) cycle
            a_soa_tend(m,ll) = (a_soa(m,ll)*tmpf - a_soa_in(m,ll))/dtfull
         end do
      end do

end subroutine modal_aero_soaexch

end module modal_aero_gasaerexch
