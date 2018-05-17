! modal_aero_gasaerexch.F90


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_gasaerexch --- does modal aerosol gas-aerosol exchange
!
! !INTERFACE:
   module modal_aero_gasaerexch

! !USES:
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use chem_mods,       only:  gas_pcnst
  use modal_aero_data, only:  nspec_max, nsoa, npoa, soa_multi_species
  use ref_pres,        only:  top_lev => clim_modal_aero_top_lev
  use ppgrid,          only:  pcols, pver
  use modal_aero_data, only:  ntot_amode, numptr_amode, sigmag_amode
  use modal_aero_data, only: lptr2_soa_g_amode, lptr2_soa_a_amode, lptr2_pom_a_amode

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public modal_aero_gasaerexch_sub, modal_aero_gasaerexch_init

! !PUBLIC DATA MEMBERS:
  integer, parameter :: pcnstxx = gas_pcnst
  integer, protected, public :: maxspec_pcage != nspec_max

  integer, protected, public :: modefrm_pcage
  integer, protected, public :: nspecfrm_pcage
  integer :: modetoo_pcage

  integer, protected, allocatable, public :: lspecfrm_pcage(:)
  integer, protected, allocatable, public :: lspectoo_pcage(:)

  real(r8), parameter, public :: n_so4_monolayers_pcage = 8.0_r8

! number of so4(+nh4) monolayers needed to "age" a carbon particle

  real(r8), parameter, public :: &
              dr_so4_monolayers_pcage = n_so4_monolayers_pcage * 4.76e-10_r8
! thickness of the so4 monolayers (m)
! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3,
!    --> 1 mol so4(+nh4)  = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
! aging criterion is approximate so do not try to distinguish
!    sulfuric acid, bisulfate, ammonium sulfate

  real(r8), protected, allocatable, public :: soa_equivso4_factor(:)
! this factor converts an soa volume to a volume of so4(+nh4)
! having same hygroscopicity as the soa

  real (r8) :: fac_m2v_nh4, fac_m2v_so4
  real (r8), allocatable :: fac_m2v_soa(:)

  real (r8), allocatable :: fac_m2v_pcarbon(:)

! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------


  contains


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
subroutine modal_aero_gasaerexch_sub(                            &
                        lchnk,    ncol,     nstep,               &
                        loffset,  deltat,                        &
                        t,        pmid,     pdel,                &
                        qh2o,               troplev,             &
                        q,                  qqcw,                &
                        dqdt_other,         dqqcwdt_other,       &
                        dgncur_a,           dgncur_awet,         &
                        sulfeq         )

! !USES:
use modal_aero_data,   only:  alnsg_amode,lmassptr_amode,cnst_name_cw
use modal_aero_data,   only:  lptr_so4_a_amode,lptr_nh4_a_amode
use modal_aero_data,   only:  modeptr_pcarbon,nspec_amode,specmw_amode,specdens_amode
use modal_aero_rename, only:  modal_aero_rename_sub

use cam_history,       only:  outfld, fieldname_len
use chem_mods,         only:  adv_mass
use constituents,      only:  pcnst, cnst_name, cnst_get_ind
use mo_tracname,       only:  solsym
use physconst,         only:  gravit, mwdry, rair
use cam_abortutils,    only:  endrun
use spmd_utils,        only:  iam, masterproc


implicit none

! !PARAMETERS:
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   integer,  intent(in)    :: troplev(pcols)       ! tropopause vertical index
   real(r8), intent(in)    :: deltat               ! time step (s)

   real(r8), intent(inout) :: q(ncol,pver,pcnstxx) ! tracer mixing ratio (TMR) array
                                                   ! *** MUST BE  #/kmol-air for number
                                                   ! *** MUST BE mol/mol-air for mass
                                                   ! *** NOTE ncol dimension
   real(r8), intent(inout) :: qqcw(ncol,pver,pcnstxx) 
                                                   ! like q but for cloud-borner tracers
   real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx) 
                                                   ! TMR tendency from other continuous
                                                   ! growth processes (aqchem, soa??)
                                                   ! *** NOTE ncol dimension
   real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx) 
                                                   ! like dqdt_other but for cloud-borner tracers
   real(r8), intent(in)    :: t(pcols,pver)        ! temperature at model levels (K)
   real(r8), intent(in)    :: pmid(pcols,pver)     ! pressure at model levels (Pa)
   real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: qh2o(pcols,pver)     ! water vapor mixing ratio (kg/kg)
   real(r8), intent(in)    :: dgncur_a(pcols,pver,ntot_amode)
   real(r8), intent(in)    :: dgncur_awet(pcols,pver,ntot_amode)
   real(r8), pointer       :: sulfeq(:,:,:)

                                 ! dry & wet geo. mean dia. (m) of number distrib.

! !DESCRIPTION: 
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
!  
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: jsrflx_gaexch = 1
   integer, parameter :: jsrflx_rename = 2
   integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
   integer, parameter :: method_soa = 2
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr modal_aero_soaexch

   integer :: i, iq, itmpa
   integer :: idiagss
   integer :: ido_so4a(ntot_amode), ido_nh4a(ntot_amode)
   integer ::  ido_soaa(ntot_amode,nsoa)
   integer :: j, jac, jsrf, jsoa
   integer :: k,p
   integer :: l, l2, lb, lsfrm, lstoo
   integer :: l_so4g, l_nh4g, l_msag
   integer :: l_soag(nsoa)
   integer :: n, nn, niter, niter_max, ntot_soamode

   logical :: is_dorename_atik, dorename_atik(ncol,pver)

   character(len=fieldname_len+3) :: fieldname
   character(len=100) :: msg !BSINGH - msg string for endrun calls

   real (r8) :: avg_uprt_nh4, avg_uprt_so4, avg_uprt_soa(nsoa)
   real (r8) :: deltatxx
   real (r8) :: dqdt_nh4(ntot_amode), dqdt_so4(ntot_amode)
   real (r8) :: dqdt_soa(ntot_amode,nsoa)
   real (r8) :: dqdt_soag(nsoa)
   real (r8) :: fac_volsfc_pcarbon
   real (r8) :: fgain_nh4(ntot_amode), fgain_so4(ntot_amode)
   real (r8) :: fgain_soa(ntot_amode,nsoa)
   real (r8) :: g0_soa(nsoa)
   real(r8)  :: mw_poa_host(npoa)          ! molec wght of poa used in host code
   real(r8)  :: mw_soa_host(nsoa)          ! molec wght of poa used in host code
   real (r8) :: pdel_fac
   real (r8) :: qmax_nh4, qnew_nh4, qnew_so4
   real (r8) :: qold_nh4(ntot_amode), qold_so4(ntot_amode)
   real (r8) :: qold_poa(ntot_amode,npoa)
   real (r8) :: qold_soa(ntot_amode,nsoa)
   real (r8) :: qold_soag(nsoa)
   real (r8) :: sum_dqdt_msa, sum_dqdt_so4
   real (r8) :: sum_dqdt_soa(nsoa)
   real (r8) :: sum_dqdt_nh4, sum_dqdt_nh4_b
   real (r8) :: sum_uprt_msa, sum_uprt_nh4, sum_uprt_so4
   real (r8) :: sum_uprt_soa(nsoa)
   real (r8) :: tmp1, tmp2, tmpa
   real (r8) :: tmp_kxt, tmp_pxt
   real (r8) :: tmp_so4a_bgn, tmp_so4a_end
   real (r8) :: tmp_so4g_avg, tmp_so4g_bgn, tmp_so4g_equ
   real (r8) :: uptkrate(ntot_amode,pcols,pver)  
   real (r8) :: uptkratebb(ntot_amode)
   real (r8) :: uptkrate_soa(ntot_amode,nsoa)  
                ! gas-to-aerosol mass transfer rates (1/s)
   real (r8) :: vol_core, vol_shell
   real (r8) :: xferfrac_pcage, xferfrac_max
   real (r8) :: xferrate

   logical  :: do_msag         ! true if msa gas is a species
   logical  :: do_nh4g         ! true if nh3 gas is a species
   logical  :: do_soag_any         ! true if soa gas is a species
   logical  :: do_soag(nsoa)       ! true if soa gas is a species

   logical  :: dotend(pcnstxx)          ! identifies species directly involved in
                                        !    gas-aerosol exchange (gas condensation)
   logical  :: dotendqqcw(pcnstxx)      ! like dotend but for cloud-borner tracers
   logical  :: dotendrn(pcnstxx), dotendqqcwrn(pcnstxx)
                                        ! identifies species involved in renaming
                                        !    after "continuous growth"
                                        !    (gas-aerosol exchange and aqchem)

   integer, parameter :: nsrflx = 2     ! last dimension of qsrflx
   real(r8) :: dqdt(ncol,pver,pcnstxx)  ! TMR "delta q" array - NOTE dims
   real(r8) :: dqqcwdt(ncol,pver,pcnstxx) ! like dqdt but for cloud-borner tracers
   real(r8) :: qsrflx(pcols,pcnstxx,nsrflx)
                              ! process-specific column tracer tendencies 
                              ! (1=renaming, 2=gas condensation)
   real(r8) :: qconff(pcols,pver),qevapff(pcols,pver)
   real(r8) :: qconbb(pcols,pver),qevapbb(pcols,pver)
   real(r8) :: qconbg(pcols,pver),qevapbg(pcols,pver)
   real(r8) :: qcon(pcols,pver),qevap(pcols,pver)

   real(r8) :: qqcwsrflx(pcols,pcnstxx,nsrflx)

!  following only needed for diagnostics
   real(r8) :: qold(ncol,pver,pcnstxx)  ! NOTE dims
   real(r8) :: qnew(ncol,pver,pcnstxx)  ! NOTE dims
   real(r8) :: qdel(ncol,pver,pcnstxx)  ! NOTE dims
   real(r8) :: dumavec(1000), dumbvec(1000), dumcvec(1000)
   real(r8) :: qqcwold(ncol,pver,pcnstxx)
   real(r8) :: dqdtsv1(ncol,pver,pcnstxx)
   real(r8) :: dqqcwdtsv1(ncol,pver,pcnstxx)


!----------------------------------------------------------------------

! set gas species indices
   call cnst_get_ind( 'H2SO4', l_so4g, .false. )
   call cnst_get_ind( 'NH3',   l_nh4g, .false. )
   call cnst_get_ind( 'MSA',   l_msag, .false. )
   l_so4g = l_so4g - loffset
   l_nh4g = l_nh4g - loffset
   l_msag = l_msag - loffset
   if ((l_so4g <= 0) .or. (l_so4g > pcnstxx)) then
      write( *, '(/a/a,2i7)' )   &
         '*** modal_aero_gasaerexch_sub -- cannot find H2SO4 species',   &
         '    l_so4g, loffset =', l_so4g, loffset
      call endrun( 'modal_aero_gasaerexch_sub error' )
   end if
   do_nh4g = .false.
   do_msag = .false.
   if ((l_nh4g > 0) .and. (l_nh4g <= pcnstxx)) do_nh4g = .true.
   if ((l_msag > 0) .and. (l_msag <= pcnstxx)) do_msag = .true.

   do_soag_any = .false.
   do_soag(:) = .false.
   do jsoa = 1, nsoa
      l_soag(jsoa) = lptr2_soa_g_amode(jsoa) - loffset
      if ((method_soa == 1) .or. (method_soa == 2)) then
         if ((l_soag(jsoa) > 0) .and. (l_soag(jsoa) <= pcnstxx)) then
            do_soag_any = .true.
            do_soag(jsoa) = .true.
         end if
      else if (method_soa /= 0) then
         write(*,'(/a,1x,i10)') '*** modal_aero_gasaerexch_sub - bad method_soa =', method_soa
         call endrun( 'modal_aero_gasaerexch_sub error' )
      end if
   end do ! jsoa

! set tendency flags
   dotend(:) = .false.
   dotendqqcw(:) = .false.
   ido_so4a(:) = 0
   ido_nh4a(:) = 0
   ido_soaa(:,:) = 0

   dotend(l_so4g) = .true.
   if ( do_nh4g ) dotend(l_nh4g) = .true.
   if ( do_msag ) dotend(l_msag) = .true.
   do jsoa = 1, nsoa
      if ( do_soag(jsoa) ) dotend(l_soag(jsoa)) = .true.
   end do

   ntot_soamode = 0
   do n = 1, ntot_amode
      l = lptr_so4_a_amode(n)-loffset
      if ((l > 0) .and. (l <= pcnstxx)) then
         dotend(l) = .true.
         ido_so4a(n) = 1
         if ( do_nh4g ) then
            l = lptr_nh4_a_amode(n)-loffset
            if ((l > 0) .and. (l <= pcnstxx)) then
               dotend(l) = .true.
               ido_nh4a(n) = 1
            end if
         end if
      end if

      do jsoa = 1, nsoa
         if ( do_soag(jsoa) ) then
            l = lptr2_soa_a_amode(n,jsoa)-loffset
            if ((l > 0) .and. (l <= pcnstxx)) then
               dotend(l) = .true.
               ido_soaa(n,jsoa) = 1
               ntot_soamode = n
            end if
         end if
      end do ! jsoa
   end do ! n


   if ( do_soag_any ) ntot_soamode = max( ntot_soamode, modefrm_pcage )

   if (modefrm_pcage > 0) then
      ido_so4a(modefrm_pcage) = 2
      if (ido_nh4a(modetoo_pcage) == 1) ido_nh4a(modefrm_pcage) = 2
      do jsoa = 1, nsoa
         if (ido_soaa(modetoo_pcage,jsoa) == 1) ido_soaa(modefrm_pcage,jsoa) = 2
      end do
      do iq = 1, nspecfrm_pcage
         lsfrm = lspecfrm_pcage(iq)-loffset
         lstoo = lspectoo_pcage(iq)-loffset
         if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
            dotend(lsfrm) = .true.
            if ((lstoo > 0) .and. (lstoo <= pcnst)) then
               dotend(lstoo) = .true.
            end if
         end if
      end do


      n = modeptr_pcarbon
      fac_volsfc_pcarbon = exp( 2.5_r8*(alnsg_amode(n)**2) )
      xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps
   end if


! zero out tendencies and other
   dqdt(:,:,:) = 0.0_r8
   dqqcwdt(:,:,:) = 0.0_r8
   qsrflx(:,:,:) = 0.0_r8
   qqcwsrflx(:,:,:) = 0.0_r8

!-------Initialize evap/cond diagnostics (ncols x pver)-----------
   qconff(:,:) = 0.0_r8
   qevapff(:,:) = 0.0_r8
   qconbb(:,:) = 0.0_r8
   qevapbb(:,:) = 0.0_r8
   qconbg(:,:) = 0.0_r8
   qevapbg(:,:) = 0.0_r8
   qcon(:,:) = 0.0_r8
   qevap(:,:) = 0.0_r8
!---------------------------------------------------

! compute gas-to-aerosol mass transfer rates
   call gas_aer_uptkrates( ncol,       loffset,                &
                           q,          t,          pmid,       &
                           dgncur_awet,            uptkrate    )


! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)


   jsrf = jsrflx_gaexch
   do k=top_lev,pver
      do i=1,ncol

!   fgain_so4(n) = fraction of total h2so4 uptake going to mode n
!   fgain_nh4(n) = fraction of total  nh3  uptake going to mode n
         sum_uprt_so4 = 0.0_r8
         sum_uprt_nh4 = 0.0_r8
         sum_uprt_soa = 0.0_r8
         do n = 1, ntot_amode
            uptkratebb(n) = uptkrate(n,i,k)
            if (ido_so4a(n) > 0) then
               fgain_so4(n) = uptkratebb(n)
               sum_uprt_so4 = sum_uprt_so4 + fgain_so4(n)
               if (ido_so4a(n) == 1) then
                  qold_so4(n) = q(i,k,lptr_so4_a_amode(n)-loffset)
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
                  qold_nh4(n) = q(i,k,lptr_nh4_a_amode(n)-loffset)
               else
                  qold_nh4(n) = 0.0_r8
               end if
            else
               fgain_nh4(n) = 0.0_r8
               qold_nh4(n) = 0.0_r8
            end if

            do j = 1, npoa
               l = lptr2_pom_a_amode(n,j)-loffset
               if (l > 0) then
                  qold_poa(n,j) = q(i,k,l)
               else
                  qold_poa(n,j) = 0.0_r8
               end if
            end do

            itmpa = 0
            do jsoa = 1, nsoa
               if (ido_soaa(n,jsoa) > 0) then
                  ! 0.81 factor is for gas diffusivity (soa/h2so4)
                  ! (differences in fuch-sutugin and accom coef ignored)
                  fgain_soa(n,jsoa) = uptkratebb(n)*0.81_r8
                  sum_uprt_soa(jsoa) = sum_uprt_soa(jsoa) + fgain_soa(n,jsoa)
                  if (ido_soaa(n,jsoa) == 1) then
                     l = lptr2_soa_a_amode(n,jsoa)-loffset
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
            do n = 1, ntot_amode
               fgain_so4(n) = fgain_so4(n) / sum_uprt_so4
            end do
         end if
!       at this point (sum_uprt_so4 <= 0.0) only when all the fgain_so4 are zero
         if (sum_uprt_nh4 > 0.0_r8) then
            do n = 1, ntot_amode
               fgain_nh4(n) = fgain_nh4(n) / sum_uprt_nh4
            end do
         end if

         do jsoa = 1, nsoa
            if (sum_uprt_soa(jsoa) > 0.0_r8) then
               do n = 1, ntot_amode
                  fgain_soa(n,jsoa) = fgain_soa(n,jsoa) / sum_uprt_soa(jsoa)
               end do
            end if
         end do

!   uptake amount (fraction of gas uptaken) over deltat
         avg_uprt_so4 = (1.0_r8 - exp(-deltatxx*sum_uprt_so4))/deltatxx
         avg_uprt_nh4 = (1.0_r8 - exp(-deltatxx*sum_uprt_nh4))/deltatxx

         do jsoa = 1, nsoa
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

         do jsoa = 1, nsoa
            if ( do_soag(jsoa) ) then
               sum_dqdt_soa(jsoa) = q(i,k,l_soag(jsoa)) * avg_uprt_soa(jsoa)
            else
               sum_dqdt_soa(jsoa) = 0.0_r8
            end if
         end do

         if ( associated(sulfeq) .and. (k <= troplev(i)) ) then
            !   compute TMR tendencies for so4 interstial aerosol due to reversible gas uptake
            !   only above the tropopause

            tmp_kxt = deltatxx*sum_uprt_so4  ! sum over modes of uptake_rate*deltat
            tmp_pxt = 0.0_r8
            do n = 1, ntot_amode
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
            do n = 1, ntot_amode
               if (ido_so4a(n) <= 0) cycle
               ! calc change to so4(a) in mode n
               if (ido_so4a(n) == 1) then
                  l = lptr_so4_a_amode(n)-loffset
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
            do n = 1, ntot_amode
               dqdt_so4(n) = fgain_so4(n)*(sum_dqdt_so4 + sum_dqdt_msa)
            end do
         end if

         !   compute TMR tendencies for nh4 interstial aerosol due to simple gas uptake
         !   but force nh4/so4 molar ratio <= 2
         sum_dqdt_nh4_b = 0.0_r8
         dqdt_nh4(:) = 0._r8
         if ( do_nh4g ) then
            do n = 1, ntot_amode
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
            do jsoa = 1, nsoa
               qold_soag(jsoa) = q(i,k,l_soag(jsoa))
            end do
            mw_poa_host = 12.0_r8
            mw_soa_host = 250.0_r8

            call modal_aero_soaexch( deltat, t(i,k), pmid(i,k), &
                 niter, niter_max, ntot_amode, ntot_soamode, npoa, nsoa, &
                 mw_poa_host, mw_soa_host, &
                 qold_soag, qold_soa, qold_poa, uptkrate_soa, &
                 dqdt_soag, dqdt_soa )
            sum_dqdt_soa(:) = -dqdt_soag(:)

         else if ( do_soag_any ) then
!   compute TMR tendencies for soa interstial aerosol
!   due to simple gas uptake

            do jsoa = 1, nsoa
               do n = 1, ntot_amode
                  dqdt_soa(n,jsoa) = fgain_soa(n,jsoa)*sum_dqdt_soa(jsoa)
               end do
            end do
         else ! method_soa is neither 1 nor 2, no uptake
            dqdt_soa(:,:) = 0.0_r8
         end if

         pdel_fac = pdel(i,k)/gravit
         do n = 1, ntot_amode
            if (ido_so4a(n) == 1) then
               l = lptr_so4_a_amode(n)-loffset
               dqdt(i,k,l) = dqdt_so4(n)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(n)*pdel_fac
            end if

            if ( do_nh4g ) then
               if (ido_nh4a(n) == 1) then
                  l = lptr_nh4_a_amode(n)-loffset
                  dqdt(i,k,l) = dqdt_nh4(n)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(n)*pdel_fac
               end if
            end if

            do jsoa = 1, nsoa
               if ( do_soag(jsoa) ) then
                  if (ido_soaa(n,jsoa) == 1) then
                     l = lptr2_soa_a_amode(n,jsoa)-loffset
                     dqdt(i,k,l) = dqdt_soa(n,jsoa) !calculated by  modal_aero_soaexch for method_soa=2
                     qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(n,jsoa)*pdel_fac
!------- Add code for condensation/evaporation diagnostics---
                     if (nsoa.eq.15) then !check for current SOA package
                        if(jsoa.ge.1.and.jsoa.le.5) then ! Fossil SOA species
                           if (dqdt_soa(n,jsoa).ge.0.0_r8) then
                              qconff(i,k)=qconff(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           elseif(dqdt_soa(n,jsoa).lt.0.0_r8) then
                              qevapff(i,k)=qevapff(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           endif

                        elseif(jsoa.ge.6.and.jsoa.le.10) then ! Biomass SOA species
                           if (dqdt_soa(n,jsoa).ge.0.0_r8) then
                              qconbb(i,k)=qconbb(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           elseif(dqdt_soa(n,jsoa).lt.0.0_r8) then
                              qevapbb(i,k)=qevapbb(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           endif

                        elseif(jsoa.ge.11.and.jsoa.le.15) then ! Biomass SOA species
                           if (dqdt_soa(n,jsoa).ge.0.0_r8) then
                              qconbg(i,k)=qconbg(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           elseif(dqdt_soa(n,jsoa).lt.0.0_r8) then
                              qevapbg(i,k)=qevapbg(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           endif

                        endif ! jsoa
                     endif !nsoa
                     if (nsoa.eq.5) then !check for current SOA package
                           if (dqdt_soa(n,jsoa).ge.0.0_r8) then
                              qcon(i,k)=qcon(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           elseif(dqdt_soa(n,jsoa).lt.0.0_r8) then
                              qevap(i,k)=qevap(i,k)+dqdt_soa(n,jsoa)*(adv_mass(l)/mwdry)
                           endif
                     endif !nsoa
!---------------------------------------------------------------------------------------------------------------------
                  end if
               end if
            end do
         end do ! n
 
!   compute TMR tendencies for h2so4, nh3, and msa gas
!   due to simple gas uptake
         l = l_so4g
         dqdt(i,k,l) = -sum_dqdt_so4
         qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac

         if ( do_msag ) then
            l = l_msag
            dqdt(i,k,l) = -sum_dqdt_msa
            qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
         end if

         if ( do_nh4g ) then
            l = l_nh4g
            dqdt(i,k,l) = -sum_dqdt_nh4_b
            qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
         end if

         do jsoa = 1, nsoa
            if ( do_soag(jsoa) ) then
               l = l_soag(jsoa)
               dqdt(i,k,l) = -sum_dqdt_soa(jsoa)
! dqdt for gas is negative of the sum of dqdt for aerosol soa species in each mode: Manish
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
            end if
         end do
 
!   compute TMR tendencies associated with primary carbon aging
         if (modefrm_pcage > 0) then
            n = modeptr_pcarbon
            tmpa = 0.0_r8
            do jsoa = 1, nsoa
               tmpa = tmpa + dqdt_soa(n,jsoa)*fac_m2v_soa(jsoa)*soa_equivso4_factor(jsoa)
            end do
            vol_shell = deltat *   &
                 ( dqdt_so4(n)*fac_m2v_so4 + dqdt_nh4(n)*fac_m2v_nh4 + tmpa )
            vol_core = 0.0_r8
            do l = 1, nspec_amode(n)
               vol_core = vol_core + &
                    q(i,k,lmassptr_amode(l,n)-loffset)*fac_m2v_pcarbon(l)
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
                  lsfrm = lspecfrm_pcage(iq)-loffset
                  lstoo = lspectoo_pcage(iq)-loffset
                  xferrate = (xferfrac_pcage/deltat)*q(i,k,lsfrm)
                  dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferrate
                  qsrflx(i,lsfrm,jsrf) = qsrflx(i,lsfrm,jsrf) - xferrate*pdel_fac
                  if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                     dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferrate
                     qsrflx(i,lstoo,jsrf) = qsrflx(i,lstoo,jsrf) + xferrate*pdel_fac
                  end if
               end do

               if (ido_so4a(modetoo_pcage) > 0) then
                  l = lptr_so4_a_amode(modetoo_pcage)-loffset
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_so4(modefrm_pcage)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(modefrm_pcage)*pdel_fac
               end if

               if (ido_nh4a(modetoo_pcage) > 0) then
                  l = lptr_nh4_a_amode(modetoo_pcage)-loffset
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_nh4(modefrm_pcage)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(modefrm_pcage)*pdel_fac
               end if

               do jsoa = 1, nsoa
                  if (ido_soaa(modetoo_pcage,jsoa) > 0) then
                     l = lptr2_soa_a_amode(modetoo_pcage,jsoa)-loffset
                     dqdt(i,k,l) = dqdt(i,k,l) + dqdt_soa(modefrm_pcage,jsoa)
                     qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(modefrm_pcage,jsoa)*pdel_fac
                  end if
               end do

            end if

         end if

      end do   ! "i = 1, ncol"
   end do     ! "k = top_lev, pver"

! set "temporary testing arrays"
   qold(:,:,:) = q(:,:,:)
   qqcwold(:,:,:) = qqcw(:,:,:)
   dqdtsv1(:,:,:) = dqdt(:,:,:)
   dqqcwdtsv1(:,:,:) = dqqcwdt(:,:,:)


!
! do renaming calcs
!
   dotendrn(:) = .false.
   dotendqqcwrn(:) = .false.
   dorename_atik(1:ncol,:) = .true.
   is_dorename_atik = .true.
   call modal_aero_rename_sub(                              &
        'modal_aero_gasaerexch_sub',            &
        lchnk,             ncol,      nstep,    &
        loffset,           deltat,              &
        pdel,              troplev,             &
        dotendrn,          q,                   &
        dqdt,              dqdt_other,          &
        dotendqqcwrn,      qqcw,                &
        dqqcwdt,           dqqcwdt_other,       &
        is_dorename_atik,  dorename_atik,       &
        jsrflx_rename,     nsrflx,              &
        qsrflx,            qqcwsrflx            )


!  This applies dqdt tendencies for all species
!  apply the dqdt to update q (and same for qqcw)
!
   do l = 1, pcnstxx
      if ( dotend(l) .or. dotendrn(l) ) then
         do k = top_lev, pver
            do i = 1, ncol
               q(i,k,l) = q(i,k,l) + dqdt(i,k,l)*deltat
            end do
         end do
      end if
      if ( dotendqqcw(l) .or. dotendqqcwrn(l) ) then
         do k = top_lev, pver
            do i = 1, ncol
               qqcw(i,k,l) = qqcw(i,k,l) + dqqcwdt(i,k,l)*deltat
            end do
         end do
      end if
   end do

! diagnostics start -------------------------------------------------------
!!$   if (ldiag3 > 0) then
!!$   if (icol_diag > 0) then 
!!$      i = icol_diag
!!$      write(*,'(a,3i5)') 'gasaerexch ppp nstep,lat,lon', nstep, latndx(i), lonndx(i)
!!$      write(*,'(2i5,3(2x,a))') 0, 0, 'ppp', 'pdel for all k'
!!$      write(*,'(1p,7e12.4)') (pdel(i,k), k=top_lev,pver)
!!$
!!$      write(*,'(a,3i5)') 'gasaerexch ddd nstep,lat,lon', nstep, latndx(i), lonndx(i)
!!$      do l = 1, pcnstxx
!!$         lb = l + loffset
!!$
!!$         if ( dotend(l) .or. dotendrn(l) ) then
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd1', cnst_name(lb),    'qold for all k'
!!$            write(*,'(1p,7e12.4)') (qold(i,k,l), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd2', cnst_name(lb),    'qnew for all k'
!!$            write(*,'(1p,7e12.4)') (q(i,k,l), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd3', cnst_name(lb),    'dqdt from conden for all k'
!!$            write(*,'(1p,7e12.4)') (dqdtsv1(i,k,l), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd4', cnst_name(lb),    'dqdt from rename for all k'
!!$            write(*,'(1p,7e12.4)') ((dqdt(i,k,l)-dqdtsv1(i,k,l)), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd5', cnst_name(lb),    'dqdt other for all k'
!!$            write(*,'(1p,7e12.4)') (dqdt_other(i,k,l), k=top_lev,pver)
!!$         end if
!!$
!!$         if ( dotendqqcw(l) .or. dotendqqcwrn(l) ) then
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd1', cnst_name_cw(lb), 'qold for all k'
!!$            write(*,'(1p,7e12.4)') (qqcwold(i,k,l), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd2', cnst_name_cw(lb), 'qnew for all k'
!!$            write(*,'(1p,7e12.4)') (qqcw(i,k,l), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd3', cnst_name_cw(lb), 'dqdt from conden for all k'
!!$            write(*,'(1p,7e12.4)') (dqqcwdtsv1(i,k,l), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd4', cnst_name_cw(lb), 'dqdt from rename for all k'
!!$            write(*,'(1p,7e12.4)') ((dqqcwdt(i,k,l)-dqqcwdtsv1(i,k,l)), k=top_lev,pver)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd5', cnst_name_cw(lb), 'dqdt other for all k'
!!$            write(*,'(1p,7e12.4)') (dqqcwdt_other(i,k,l), k=top_lev,pver)
!!$         end if
!!$
!!$      end do
!!$
!!$      write(*,'(a,3i5)') 'gasaerexch fff nstep,lat,lon', nstep, latndx(i), lonndx(i)
!!$      do l = 1, pcnstxx
!!$         lb = l + loffset
!!$         if ( dotend(l) .or. dotendrn(l) .or. dotendqqcw(l) .or. dotendqqcwrn(l) ) then
!!$            write(*,'(i5,2(2x,a,2l3))') l, &
!!$               cnst_name(lb), dotend(l), dotendrn(l), &
!!$               cnst_name_cw(lb), dotendqqcw(l), dotendqqcwrn(l)
!!$         end if
!!$      end do
!!$
!!$   end if
!!$   end if
! diagnostics end ---------------------------------------------------------

!-----Outfld for condensation/evaporation------------------------------
   if (nsoa.eq.5) then !check for current SOA package
      call outfld(trim('qcon_gaex'), qcon(:,:), pcols, lchnk )
      call outfld(trim('qevap_gaex'), qevap(:,:), pcols, lchnk )
   endif
!-----------------------------------------------------------------------
   if (nsoa.eq.15) then !check for current SOA package
      call outfld(trim('qconff_gaex'), qconff(:,:), pcols, lchnk )
      call outfld(trim('qevapff_gaex'), qevapff(:,:), pcols, lchnk )
      call outfld(trim('qconbb_gaex'), qconbb(:,:), pcols, lchnk )
      call outfld(trim('qevapbb_gaex'), qevapbb(:,:), pcols, lchnk )
      call outfld(trim('qconbg_gaex'), qconbg(:,:), pcols, lchnk )
      call outfld(trim('qevapbg_gaex'), qevapbg(:,:), pcols, lchnk )
   endif
!-----------------------------------------------------------------------

!   do history file column-tendency fields
   do l = 1, pcnstxx
      lb = l + loffset
      do jsrf = 1, 2
         do jac = 1, 2
	    if (jac == 1) then
               if (jsrf == jsrflx_gaexch) then
                  if ( .not. dotend(l) ) cycle
                  fieldname = trim(cnst_name(lb)) // '_sfgaex1'
               else if (jsrf == jsrflx_rename) then
                  if ( .not. dotendrn(l) ) cycle
                  fieldname = trim(cnst_name(lb)) // '_sfgaex2'
               else 
                  cycle
               end if
               do i = 1, ncol
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do
               call outfld( fieldname, qsrflx(:,l,jsrf), pcols, lchnk )
	    else
               if (jsrf == jsrflx_gaexch) then
                  cycle
               else if (jsrf == jsrflx_rename) then
                  if ( .not. dotendqqcwrn(l) ) cycle
                  fieldname = trim(cnst_name_cw(lb)) // '_sfgaex2'
               else 
                  cycle
               end if
               do i = 1, ncol
                  qqcwsrflx(i,l,jsrf) = qqcwsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do
               call outfld( fieldname, qqcwsrflx(:,l,jsrf), pcols, lchnk )
	    end if
         end do ! jac = ...
      end do ! jsrf = ...
   end do ! l = ...
   
   return
   end subroutine modal_aero_gasaerexch_sub


!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine gas_aer_uptkrates( ncol,       loffset,                &
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

use physconst, only: mwdry, rair

implicit none


   integer,  intent(in) :: ncol                 ! number of atmospheric column
   integer,  intent(in) :: loffset
   real(r8), intent(in) :: q(ncol,pver,pcnstxx) ! Tracer array (mol,#/mol-air)
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature in Kelvin
   real(r8), intent(in) :: pmid(pcols,pver)     ! Air pressure in Pa
   real(r8), intent(in) :: dgncur_awet(pcols,pver,ntot_amode)

   real(r8), intent(out) :: uptkrate(ntot_amode,pcols,pver)  
                            ! gas-to-aerosol mass transfer rates (1/s)


! local 
   integer, parameter :: nghq = 2
   integer :: i, iq, k, l1, l2, la, n

   ! Can use sqrt here once Lahey is gone.
   real(r8), parameter :: tworootpi = 3.5449077_r8
   real(r8), parameter :: root2 = 1.4142135_r8
   real(r8), parameter :: beta = 2.0_r8

   real(r8) :: aircon
   real(r8) :: const
   real(r8) :: dp, dum_m2v
   real(r8) :: dryvol_a(pcols,pver)
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
   do n = 1, ntot_amode

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

         rhoair = pmid(i,k)/(rair*t(i,k))   ! (kg-air/m3)
!        aircon = 1.0e3*rhoair/mwdry        ! (mol-air/m3)

!!   "bounded" number conc. (#/m3)
!        num_a = dryvol_a(i,k)*v2ncur_a(i,k,n)*aircon

!   number conc. (#/m3) -- note q(i,k,numptr) is (#/kmol-air) 
!   so need aircon in (kmol-air/m3)
         aircon = rhoair/mwdry              ! (kmol-air/m3)
         num_a = q(i,k,numptr_amode(n)-loffset)*aircon

!   gasdiffus = h2so4 gas diffusivity from mosaic code (m^2/s)
!               (pmid must be Pa)
         gasdiffus = 0.557e-4_r8 * (t(i,k)**1.75_r8) / pmid(i,k)
!   gasspeed = h2so4 gas mean molecular speed from mosaic code (m/s)
         gasspeed  = 1.470e1_r8 * sqrt(t(i,k))
!   freepathx2 = 2 * (h2so4 mean free path)  (m)
         freepathx2 = 6.0_r8*gasdiffus/gasspeed

         lnsg   = log( sigmag_amode(n) )
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


   return
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

      use mo_constants, only: rgas ! Gas constant (J/K/mol)

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

! Changed by Manish Shrivastava
      opoa_frac(:) = 0.0_r8 !POA does not form solution with SOA for all runs; set opoa_frac=0.0_r8  by Manish Shrivastava
      mw_poa(:) = 250.0_r8
      mw_soa(:) = 250.0_r8

      ! New SOA properties added by Manish Shrivastava on 09/27/2012
      if (ntot_soaspec ==1) then
         p0_soa_298(:) = 1.0e-10_r8
         delh_vap_soa(:) = 156.0e3_r8
         opoa_frac(:) = 0.1_r8
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

      tmpf = mw_poa_host(1)/mw_poa(1)
      do m = 1, ntot_soamode
         if ( skip_soamode(m) ) cycle
         a_opoa(m) = 0.0_r8
         do ll = 1, ntot_poaspec
            tmpf = mw_poa_host(ll)/mw_poa(ll)
            a_opoa(m) = opoa_frac(ll)*a_poa_in(m,ll)
            a_opoa(m) = max( a_opoa(m), 1.0e-20_r8 )  ! force to small non-zero value
         end do
      end do

      ! calc ambient equilibrium soa gas
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
              exp( -(delh_vap_soa(ll)/rgas)*((1.0_r8/temp)-(1.0_r8/298.0_r8)) )
         g0_soa(ll) = 1.01325e5_r8*p0_soa(ll)/pres
      end do
      ! IF mw of soa EQ 12 (as in the MAM3 default case), this has to be in
      ! should actully talk the mw from the chemistry mechanism and substitute with 12.0
      if (.not.soa_multi_species) then
         g0_soa = g0_soa*(150.0_r8/12.0_r8)
      else
      end if

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


      return

      end subroutine modal_aero_soaexch

!----------------------------------------------------------------------

!----------------------------------------------------------------------

      subroutine modal_aero_gasaerexch_init

!-----------------------------------------------------------------------
!
! Purpose:
!    set do_adjust and do_aitken flags
!    create history fields for column tendencies associated with
!       modal_aero_calcsize
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

use modal_aero_data
use modal_aero_rename

use cam_abortutils, only: endrun
use cam_history,    only: addfld, add_default, fieldname_len, horiz_only
use constituents,   only: pcnst, cnst_get_ind, cnst_name
use spmd_utils,     only: masterproc
use phys_control,   only: phys_getopts

implicit none

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
   integer  :: ipair, iq, iqfrm, iqfrm_aa, iqtoo, iqtoo_aa
   integer  :: jac,jsoa,p
   integer  :: l, l1, l2, lsfrm, lstoo, lunout
   integer  :: l_so4g, l_nh4g, l_msag
   integer  :: m, mfrm, mtoo
   integer  :: n, nacc, nait
   integer  :: nchfrm, nchfrmskip, nchtoo, nchtooskip, nspec

   logical  :: do_msag, do_nh4g
   logical  :: do_soag_any, do_soag(nsoa)
   logical  :: dotend(pcnst), dotendqqcw(pcnst)

   real(r8) :: tmp1, tmp2

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(128)                 :: msg
   character(8)                   :: unit

   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
   logical                        :: history_aerocom    ! Output the aerocom history
   !-----------------------------------------------------------------------
 
        call phys_getopts( history_aerosol_out        = history_aerosol   )
       
        maxspec_pcage = nspec_max
        allocate(lspecfrm_pcage(maxspec_pcage))
        allocate(lspectoo_pcage(maxspec_pcage))
        allocate(soa_equivso4_factor(nsoa))
        allocate(fac_m2v_soa(nsoa))
        allocate(fac_m2v_pcarbon(nspec_max))
	lunout = 6
!
!   define "from mode" and "to mode" for primary carbon aging
!
!   skip (turn off) aging if either is absent, 
!      or if accum mode so4 is absent
!
	modefrm_pcage = -999888777
	modetoo_pcage = -999888777
	if ((modeptr_pcarbon <= 0) .or. (modeptr_accum <= 0)) goto 15000
	l = lptr_so4_a_amode(modeptr_accum)
	if ((l < 1) .or. (l > pcnst)) goto 15000

	modefrm_pcage = modeptr_pcarbon
	modetoo_pcage = modeptr_accum

!
!   define species involved in each primary carbon aging pairing
!	(include aerosol water)
!   
!
	mfrm = modefrm_pcage
	mtoo = modetoo_pcage

	if (mfrm < 10) then
	    nchfrmskip = 1
	else if (mfrm < 100) then
	    nchfrmskip = 2
	else
	    nchfrmskip = 3
	end if
	if (mtoo < 10) then
	    nchtooskip = 1
	else if (mtoo < 100) then
	    nchtooskip = 2
	else
	    nchtooskip = 3
	end if
	nspec = 0

aa_iqfrm: do iqfrm = -1, nspec_amode(mfrm)

	    if (iqfrm == -1) then
		lsfrm = numptr_amode(mfrm)
		lstoo = numptr_amode(mtoo)
	    else if (iqfrm == 0) then
!   bypass transfer of aerosol water due to primary-carbon aging
		cycle aa_iqfrm
!               lsfrm = lwaterptr_amode(mfrm)
!               lstoo = lwaterptr_amode(mtoo)
	    else
		lsfrm = lmassptr_amode(iqfrm,mfrm)
		lstoo = 0
	    end if
	    if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle aa_iqfrm

	    if (lsfrm>0 .and. iqfrm>0 ) then
		nchfrm = len( trim( cnst_name(lsfrm) ) ) - nchfrmskip

! find "too" species having same lspectype_amode as the "frm" species
! AND same cnst_name (except for last 1/2/3 characters which are the mode index)
		do iqtoo = 1, nspec_amode(mtoo)
!		    if ( lspectype_amode(iqtoo,mtoo) .eq.   &
!		         lspectype_amode(iqfrm,mfrm) ) then
			lstoo = lmassptr_amode(iqtoo,mtoo)
			nchtoo = len( trim( cnst_name(lstoo) ) ) - nchtooskip
			if (cnst_name(lsfrm)(1:nchfrm) == cnst_name(lstoo)(1:nchtoo)) then
			    exit
			else
			    lstoo = 0
			end if
!		    end if
		end do
	    end if

	    if ((lstoo < 1) .or. (lstoo > pcnst)) lstoo = 0
	    nspec = nspec + 1
	    lspecfrm_pcage(nspec) = lsfrm
	    lspectoo_pcage(nspec) = lstoo
	end do aa_iqfrm

	nspecfrm_pcage = nspec

!
!   output results
!
	if ( masterproc ) then

	write(lunout,9310)

	  mfrm = modefrm_pcage
	  mtoo = modetoo_pcage
	  write(lunout,9320) 1, mfrm, mtoo

	  do iq = 1, nspecfrm_pcage
	    lsfrm = lspecfrm_pcage(iq)
	    lstoo = lspectoo_pcage(iq)
	    if (lstoo .gt. 0) then
		write(lunout,9330) lsfrm, cnst_name(lsfrm),   &
      			lstoo, cnst_name(lstoo)
	    else
		write(lunout,9340) lsfrm, cnst_name(lsfrm)
	    end if
	  end do

	write(lunout,*)

	end if ! ( masterproc ) 

9310	format( / 'subr. modal_aero_gasaerexch_init - primary carbon aging pointers' )
9320	format( 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )


15000 continue

! set tendency flags and gas species indices and flags
      dotend(:) = .false.

      call cnst_get_ind( 'H2SO4', l_so4g, .false. )
      if ((l_so4g <= 0) .or. (l_so4g > pcnst)) then
         write( *, '(/a/a,2i7)' )   &
            '*** modal_aero_gasaerexch_init -- cannot find H2SO4 species',   &
            '    l_so4g=', l_so4g
         call endrun( 'modal_aero_gasaerexch_init error' )
      end if
      dotend(l_so4g) = .true.

      call cnst_get_ind( 'NH3',   l_nh4g, .false. )
      do_nh4g = .false.
      if ((l_nh4g > 0) .and. (l_nh4g <= pcnst)) then
         do_nh4g = .true.
         dotend(l_nh4g) = .true.
      end if

      call cnst_get_ind( 'MSA',   l_msag, .false. )
      do_msag = .false.
      if ((l_msag > 0) .and. (l_msag <= pcnst)) then
         do_msag = .true.
         dotend(l_msag) = .true.
      end if

      do_soag_any = .false.
      do_soag(:) = .false.
      do jsoa = 1, nsoa
         l = lptr2_soa_g_amode(jsoa)
         if ((l > 0) .and. (l <= pcnst)) then
            do_soag_any = .true.
            do_soag(jsoa) = .true.
            dotend(l) = .true.
         end if
      end do


      do n = 1, ntot_amode
         l = lptr_so4_a_amode(n)
         if ((l > 0) .and. (l <= pcnst)) then
            dotend(l) = .true.
            if ( do_nh4g ) then
               l = lptr_nh4_a_amode(n)
               if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
            end if
         end if
         do jsoa = 1, nsoa
            if ( do_soag(jsoa) ) then
               l = lptr2_soa_a_amode(n,jsoa)
               if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
            end if
         end do
      end do

      if (modefrm_pcage > 0) then
         do iq = 1, nspecfrm_pcage
            lsfrm = lspecfrm_pcage(iq)
            lstoo = lspectoo_pcage(iq)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotend(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotend(lstoo) = .true.
               end if
            end if
         end do
      end if

!---------define history fields for new cond/evap diagnostics----------------------------------------
      fieldname=trim('qconff_gaex')
      long_name = trim('3D fields for Fossil SOA condensation')
      unit = 'kg/kg/s'
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qconff addfld', fieldname, unit

      fieldname=trim('qevapff_gaex')
      long_name = trim('3D fields for Fossil SOA evaporation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qevapff addfld', fieldname, unit

      fieldname=trim('qconbb_gaex')
      long_name = trim('3D fields for Biomass SOA condensation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qconbb addfld', fieldname, unit

      fieldname=trim('qevapbb_gaex')
      long_name = trim('3D fields for Biomass SOA evaporation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qevapbb addfld', fieldname, unit

      fieldname=trim('qconbg_gaex')
      long_name = trim('3D fields for Biogenic SOA condensation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qconbg addfld', fieldname, unit

      fieldname=trim('qevapbg_gaex')
      long_name = trim('3D fields for Biogenic SOA evaporation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qevapbg addfld', fieldname, unit 

      fieldname=trim('qcon_gaex')
      long_name = trim('3D fields for SOA condensation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qcon addfld', fieldname, unit

      fieldname=trim('qevap_gaex')
      long_name = trim('3D fields for Biogenic SOA evaporation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
      if ( masterproc ) write(*,'(3(a,3x))') 'qevap addfld', fieldname, unit 
!------------------------------------------------------------------------------

!  define history fields for basic gas-aer exchange
!  and primary carbon aging from that
      do l = 1, pcnst
         if ( .not. dotend(l) ) cycle

         tmpnamea = cnst_name(l)
         fieldname = trim(tmpnamea) // '_sfgaex1'
         long_name = trim(tmpnamea) // ' gas-aerosol-exchange primary column tendency'
         unit = 'kg/m2/s'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if ( history_aerosol ) then 
            call add_default( fieldname, 1, ' ' )
         endif
         if ( masterproc ) write(*,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit

      end do   ! l = ...
!  define history fields for aitken-->accum renaming
      dotend(:) = .false.
      dotendqqcw(:) = .false.
      do ipair = 1, npair_renamexf
         do iq = 1, nspecfrm_renamexf(ipair)
            lsfrm = lspecfrma_renamexf(iq,ipair)
            lstoo = lspectooa_renamexf(iq,ipair)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotend(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotend(lstoo) = .true.
               end if
            end if

            lsfrm = lspecfrmc_renamexf(iq,ipair)
            lstoo = lspectooc_renamexf(iq,ipair)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotendqqcw(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotendqqcw(lstoo) = .true.
               end if
            end if
         end do ! iq = ...
      end do ! ipair = ...

      do l = 1, pcnst
      do jac = 1, 2
         if (jac == 1) then
            if ( .not. dotend(l) ) cycle
            tmpnamea = cnst_name(l)
         else
            if ( .not. dotendqqcw(l) ) cycle
            tmpnamea = cnst_name_cw(l)
         end if

         fieldname = trim(tmpnamea) // '_sfgaex2'
         long_name = trim(tmpnamea) // ' gas-aerosol-exchange renaming column tendency'
         unit = 'kg/m2/s'
         if ((tmpnamea(1:3) == 'num') .or. &
             (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if ( history_aerosol ) then 
            call add_default( fieldname, 1, ' ' )
         endif
         if ( masterproc ) write(*,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
      end do   ! jac = ...
      end do   ! l = ...


! set for used in aging calcs:  
!    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
!    soa_equivso4_factor(:)
      soa_equivso4_factor = 0.0_r8
      if (modefrm_pcage > 0) then
         n = modeptr_accum
         l = lptr_so4_a_amode(n) ; l2 = -1
         if (l <= 0) call endrun( 'modal_aero_gasaerexch_init error a001 finding accum. so4' )
         do l1 = 1, nspec_amode(n)
            if (lmassptr_amode(l1,n) == l) then
!               l2 = lspectype_amode(l1,n)
               l2 = l1
!               fac_m2v_so4 = specmw_amode(l2) / specdens_amode(l2)
               fac_m2v_so4 = specmw_amode(l1,n) / specdens_amode(l1,n)
!               tmp2 = spechygro(l2)
               tmp2 = spechygro(l1,n)

            end if
         end do
         if (l2 <= 0) call endrun( 'modal_aero_gasaerexch_init error a002 finding accum. so4' )

         l = lptr_nh4_a_amode(n) ; l2 = -1
         if (l > 0) then
            do l1 = 1, nspec_amode(n)
               if (lmassptr_amode(l1,n) == l) then
!                  l2 = lspectype_amode(l1,n)
                  l2 = l1
!                  fac_m2v_nh4 = specmw_amode(l2) / specdens_amode(l2)
                  fac_m2v_nh4 = specmw_amode(l1,n) / specdens_amode(l1,n)

               end if
            end do
            if (l2 <= 0) call endrun( 'modal_aero_gasaerexch_init error a002 finding accum. nh4' )
         else
            fac_m2v_nh4 = fac_m2v_so4
         end if

         do jsoa = 1, nsoa
            l = lptr2_soa_a_amode(n,jsoa) ; l2 = -1
            if (l <= 0) then
               write( msg, '(a,i4)') 'modal_aero_gasaerexch_init error a001 finding accum. jsoa =', jsoa
               call endrun( msg )
            end if
            do l1 = 1, nspec_amode(n)
               if (lmassptr_amode(l1,n) == l) then
!                  l2 = lspectype_amode(l1,n)
                  l2 = l1
!                  fac_m2v_soa(jsoa) = specmw_amode(l2) / specdens_amode(l2)
                  fac_m2v_soa(jsoa) = specmw_amode(l1,n) / specdens_amode(l1,n)
!                  soa_equivso4_factor(jsoa) = spechygro(l2)/tmp2
                  soa_equivso4_factor(jsoa) = spechygro(l1,n)/tmp2
               end if
            end do
            if (l2 <= 0) then
               write( msg, '(a,i4)') 'modal_aero_gasaerexch_init error a002 finding accum. jsoa =', jsoa
               call endrun( msg )
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


      return

      end subroutine modal_aero_gasaerexch_init


!----------------------------------------------------------------------

end module modal_aero_gasaerexch

