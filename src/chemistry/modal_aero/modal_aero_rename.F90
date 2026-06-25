! Portable code for modal aerosol mode merging (renaming)
module modal_aero_rename
  use shr_kind_mod,    only: r8 => shr_kind_r8

  implicit none
  private
  save

  public :: modal_aero_rename_init
  public :: modal_aero_rename_run

  integer, parameter, public :: maxpair_renamexf = 3
  integer, parameter, public :: method_optbb_renamexf = 2

  ! rename's OWN precomputed physics coefficients (accum-coarse-exchange path) are
  ! set by modal_aero_rename_init and read by modal_aero_rename_run.

  integer :: iulog
  character(len=32), allocatable :: cnst_name(:), cnst_name_cw(:)

  ! precomputed physics coefficients (accum-coarse-exchange path)
  integer, allocatable :: ido_mode_calcaa(:)
  real (r8) :: dp_belowcut(maxpair_renamexf)
  real (r8) :: dp_cut(maxpair_renamexf)
  real (r8) :: dp_xferall_thresh(maxpair_renamexf)
  real (r8) :: dp_xfernone_threshaa(maxpair_renamexf)
  real (r8), allocatable :: dryvol_smallest(:)
  real (r8), allocatable :: factoraa(:)
  real (r8), allocatable :: factoryy(:)
  real (r8) :: lndp_cut(maxpair_renamexf)
  real (r8) :: factor_3alnsg2(maxpair_renamexf)
  real (r8), allocatable :: v2nhirlx(:), v2nlorlx(:)
contains

  !------------------------------------------------------------------
  ! Precompute rename's own accum-coarse-exchange physics coefficients from the
  ! host-provided mode metadata + resolved pair tables (passed as arguments),
  ! and store cnst_name/iulog for diagnostics.  The shared tables/metadata are
  ! NOT retained; modal_aero_rename_run receives them as arguments each call.
  !------------------------------------------------------------------
  subroutine modal_aero_rename_init(                             &
       modal_accum_coarse_exch,                                  &
       ntot_amode,                                               &
       alnsg_amode,        dgnum_amode,                          &
       dgnumhi_amode,      dgnumlo_amode,                        &
       voltonumblo_amode,  voltonumbhi_amode,                    &
       modeptr_accum,      modeptr_coarse,    modeptr_stracoar,  &
       npair_renamexf,     modefrm_renamexf,  modetoo_renamexf,  &
       nspecfrm_renamexf,                                        &
       lspecfrma_renamexf, lspecfrmc_renamexf,                   &
       lspectooa_renamexf, lspectooc_renamexf,                   &
       igrow_shrink_renamexf, ixferable_all_renamexf,            &
       cnst_name_in,       cnst_name_cw_in,                      &
       pi,                 amRoot,            iulog_in,          &
       errmsg,             errflg                                )

    ! arguments
    logical,  intent(in) :: modal_accum_coarse_exch
    integer,  intent(in) :: ntot_amode
    real(r8), intent(in) :: alnsg_amode(:)
    real(r8), intent(in) :: dgnum_amode(:)
    real(r8), intent(in) :: dgnumhi_amode(:)
    real(r8), intent(in) :: dgnumlo_amode(:)
    real(r8), intent(in) :: voltonumblo_amode(:)
    real(r8), intent(in) :: voltonumbhi_amode(:)
    integer,  intent(in) :: modeptr_accum
    integer,  intent(in) :: modeptr_coarse
    integer,  intent(in) :: modeptr_stracoar
    integer,  intent(in) :: npair_renamexf
    integer,  intent(in) :: modefrm_renamexf(:)
    integer,  intent(in) :: modetoo_renamexf(:)
    integer,  intent(in) :: nspecfrm_renamexf(:)
    integer,  intent(in) :: lspecfrma_renamexf(:,:)
    integer,  intent(in) :: lspecfrmc_renamexf(:,:)
    integer,  intent(in) :: lspectooa_renamexf(:,:)
    integer,  intent(in) :: lspectooc_renamexf(:,:)
    integer,  intent(in) :: igrow_shrink_renamexf(:)
    integer,  intent(in) :: ixferable_all_renamexf(:)
    character(len=*), intent(in) :: cnst_name_in(:)
    character(len=*), intent(in) :: cnst_name_cw_in(:)
    real(r8), intent(in) :: pi
    logical,  intent(in) :: amRoot
    integer,  intent(in) :: iulog_in
    character(len=*), intent(out) :: errmsg
    integer,  intent(out) :: errflg

    ! local (used by the precompute + one-time log below)
    integer  :: ipair, iq, lsfrma, lsfrmc, lstooa, lstooc, lunout
    integer  :: mfrm, mtoo
    real(r8) :: frelax
    logical  :: masterproc

    errmsg = ''
    errflg = 0

    iulog = iulog_in
    ! stored only for the (disabled) per-column diagnostic + one-time log
    allocate(cnst_name(size(cnst_name_in)));       cnst_name(:)    = cnst_name_in(:)
    allocate(cnst_name_cw(size(cnst_name_cw_in))); cnst_name_cw(:) = cnst_name_cw_in(:)

    allocate(ido_mode_calcaa(ntot_amode))
    allocate(dryvol_smallest(ntot_amode))
    allocate(factoraa(ntot_amode))
    allocate(factoryy(ntot_amode))
    allocate(v2nhirlx(ntot_amode), v2nlorlx(ntot_amode))

    ! nothing to precompute unless there are renaming pairs
    if (npair_renamexf .le. 0) return

    lunout    = iulog
    masterproc = amRoot

    if (modal_accum_coarse_exch) then

!
!
!   initialize some working variables
!
!
  ido_mode_calcaa(:) = 0
  frelax = 27.0_r8

  do ipair = 1, npair_renamexf
      mfrm = modefrm_renamexf(ipair)
      mtoo = modetoo_renamexf(ipair)
      ido_mode_calcaa(mfrm) = 1

      factoraa(mfrm) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mfrm)**2))
      factoraa(mtoo) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mtoo)**2))
      factoryy(mfrm) = sqrt( 0.5_r8 )/alnsg_amode(mfrm)

!   dryvol_smallest is a very small volume mixing ratio (m3-AP/kmol-air)
!   used for avoiding overflow.  it corresponds to dp = 1 nm
!   and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air
            dryvol_smallest(mfrm) = 1.0e-25_r8
            v2nlorlx(mfrm) = voltonumblo_amode(mfrm)*frelax
            v2nhirlx(mfrm) = voltonumbhi_amode(mfrm)/frelax

      factor_3alnsg2(ipair) = 3.0_r8 * (alnsg_amode(mfrm)**2)

      dp_cut(ipair) = sqrt( &
                 dgnum_amode(mfrm)*exp(1.5_r8*(alnsg_amode(mfrm)**2)) *   &
                 dgnum_amode(mtoo)*exp(1.5_r8*(alnsg_amode(mtoo)**2)) )
      dp_xferall_thresh(ipair) = dgnum_amode(mtoo)
      dp_xfernone_threshaa(ipair) = dgnum_amode(mfrm)
      if (((mfrm == modeptr_accum) .and. (mtoo == modeptr_coarse)).or.&
                ((mfrm == modeptr_accum) .and. (mtoo == modeptr_stracoar))) then
               dp_cut(ipair)               = 4.4e-7_r8
               dp_xfernone_threshaa(ipair) = 1.6e-7_r8
               dp_xferall_thresh(ipair)    = 4.7e-7_r8
            else if (((mfrm == modeptr_coarse) .and. (mtoo == modeptr_accum)).or.&
                     ((mfrm == modeptr_stracoar) .and. (mtoo == modeptr_accum))) then
               dp_cut(ipair)               = 4.4e-7_r8
               dp_xfernone_threshaa(ipair) = 4.4e-7_r8
               dp_xferall_thresh(ipair)    = 4.1e-7_r8
            end if

      lndp_cut(ipair) = log( dp_cut(ipair) )
      dp_belowcut(ipair) = 0.99_r8*dp_cut(ipair)
         end do

!
!   output results
!
  if ( masterproc ) then

  write(lunout,9310)
  write(lunout,'(a,1x,i12)') 'method_optbb_renamexf', method_optbb_renamexf

  do 2900 ipair = 1, npair_renamexf
  mfrm = modefrm_renamexf(ipair)
  mtoo = modetoo_renamexf(ipair)
  write(lunout,9320) ipair, mfrm, mtoo, &
      igrow_shrink_renamexf(ipair), ixferable_all_renamexf(ipair)

  do iq = 1, nspecfrm_renamexf(ipair)
      lsfrma = lspecfrma_renamexf(iq,ipair)
      lstooa = lspectooa_renamexf(iq,ipair)
      lsfrmc = lspecfrmc_renamexf(iq,ipair)
      lstooc = lspectooc_renamexf(iq,ipair)
      if (lstooa .gt. 0) then
    write(lunout,9330) lsfrma, cnst_name(lsfrma),   &
           lstooa, cnst_name(lstooa)
      else
    write(lunout,9340) lsfrma, cnst_name(lsfrma)
      end if
      if (lstooc .gt. 0) then
    write(lunout,9330) lsfrmc, cnst_name_cw(lsfrmc),   &
           lstooc, cnst_name_cw(lstooc)
      else if (lsfrmc .gt. 0) then
    write(lunout,9340) lsfrmc, cnst_name_cw(lsfrmc)
      else
    write(lunout,9350)
      end if
  end do

  if (igrow_shrink_renamexf(ipair) > 0) then
    write(lunout,'(5x,a,1p,2e12.3)') 'mfrm dgnum, dgnumhi ', &
      dgnum_amode(mfrm), dgnumhi_amode(mfrm)
    write(lunout,'(5x,a,1p,2e12.3)') 'mtoo dgnum, dgnumlo ', &
      dgnum_amode(mtoo), dgnumlo_amode(mtoo)
  else
    write(lunout,'(5x,a,1p,2e12.3)') 'mfrm dgnum, dgnumlo ', &
      dgnum_amode(mfrm), dgnumlo_amode(mfrm)
    write(lunout,'(5x,a,1p,2e12.3)') 'mtoo dgnum, dgnumhi ', &
      dgnum_amode(mtoo), dgnumhi_amode(mtoo)
  end if

  write(lunout,'(5x,a,1p,2e12.3)') 'dp_cut              ', &
    dp_cut(ipair)
  write(lunout,'(5x,a,1p,2e12.3)') 'dp_xfernone_threshaa', &
    dp_xfernone_threshaa(ipair)
  write(lunout,'(5x,a,1p,2e12.3)') 'dp_xferall_thresh   ', &
    dp_xferall_thresh(ipair)

2900  continue
  write(lunout,*)

  end if ! ( masterproc )
    end if

    return

9310  format( / 'subr. modal_aero_rename_acc_crs_init' )
9320  format( / 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3, &
          5x, 'igrow_shrink', i3, 5x, 'ixferable_all', i3 )
9330  format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340  format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350  format( 5x, 'no corresponding activated species' )

  end subroutine modal_aero_rename_init

  subroutine modal_aero_rename_run(                      &
       ncol,                                   &
       loffset,           deltat,              &
       pdel,              troplev,             &
       dotendrn,          q,                   &
       dqdt,              dqdt_other,          &
       dotendqqcwrn,      qqcw,                &
       dqqcwdt,           dqqcwdt_other,       &
       is_dorename_atik,  dorename_atik,       &
       jsrflx_rename,     nsrflx,              &
       qsrflx,            qqcwsrflx,           &
       dqdt_rnpos,                             &
       ntot_amode,        npair_renamexf,      &
       modefrm_renamexf,  modetoo_renamexf,    &
       nspecfrm_renamexf,                      &
       lspecfrma_renamexf, lspecfrmc_renamexf, &
       lspectooa_renamexf, lspectooc_renamexf, &
       alnsg_amode,       voltonumblo_amode,   &
       voltonumbhi_amode, dgnum_amode,         &
       nspec_amode,       specmw_amode,        &
       specdens_amode,    lmassptr_amode,      &
       lmassptrcw_amode,  numptr_amode,        &
       numptrcw_amode,    pi,                  &
       modeptr_accum,     modeptr_coarse,      &
       modeptr_stracoar,                       &
       igrow_shrink_renamexf,                  &
       ixferable_all_renamexf,                 &
       ixferable_a_renamexf, ixferable_c_renamexf, &
       strat_only_renamexf,                    &
       modal_accum_coarse_exch,                &
       pver,              gravit,              &
       errmsg,            errflg               )
    integer,  intent(in)    :: ncol                 ! number of atmospheric column
    integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
    real(r8), intent(in)    :: deltat               ! time step (s)
    integer,  intent(in)    :: troplev(:)
    real(r8), intent(in)    :: pdel(:,:)            ! pressure thickness of levels (Pa)
    real(r8), intent(in)    :: q(:,:,:)             ! tracer mixing ratio array (mol/mol-air or #/mol-air)
    real(r8), intent(in)    :: qqcw(:,:,:)          ! like q but for cloud-borne species
    real(r8), intent(inout) :: dqdt(:,:,:)          ! TMR tendency array (renaming tendencies added on)
    real(r8), intent(inout) :: dqqcwdt(:,:,:)
    real(r8), intent(in)    :: dqdt_other(:,:,:)    ! tendencies for "other" continuous growth process
    real(r8), intent(in)    :: dqqcwdt_other(:,:,:)
    logical,  intent(inout) :: dotendrn(:)          ! species with renaming dqdt computed
    logical,  intent(inout) :: dotendqqcwrn(:)
    logical,  intent(in)    :: is_dorename_atik     ! true if dorename_atik is provided
    logical,  intent(in)    :: dorename_atik(:,:)   ! true if renaming should be done at i,k
    integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
    integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx
    real(r8), intent(out)   :: qsrflx(:,:,:)        ! process-specific column tracer tendencies
    real(r8), intent(out)   :: qqcwsrflx(:,:,:)
    real(r8), intent(out)   :: dqdt_rnpos(:,:,:)  ! positive (production) part of renaming tendency
    ! shared mode metadata + resolved renaming-pair tables (host-owned; passed in)
    integer,  intent(in)    :: ntot_amode
    integer,  intent(in)    :: npair_renamexf
    integer,  intent(in)    :: modefrm_renamexf(:)
    integer,  intent(in)    :: modetoo_renamexf(:)
    integer,  intent(in)    :: nspecfrm_renamexf(:)
    integer,  intent(in)    :: lspecfrma_renamexf(:,:)
    integer,  intent(in)    :: lspecfrmc_renamexf(:,:)
    integer,  intent(in)    :: lspectooa_renamexf(:,:)
    integer,  intent(in)    :: lspectooc_renamexf(:,:)
    real(r8), intent(in)    :: alnsg_amode(:)
    real(r8), intent(in)    :: voltonumblo_amode(:)
    real(r8), intent(in)    :: voltonumbhi_amode(:)
    real(r8), intent(in)    :: dgnum_amode(:)
    integer,  intent(in)    :: nspec_amode(:)
    real(r8), intent(in)    :: specmw_amode(:,:)
    real(r8), intent(in)    :: specdens_amode(:,:)
    integer,  intent(in)    :: lmassptr_amode(:,:)
    integer,  intent(in)    :: lmassptrcw_amode(:,:)
    integer,  intent(in)    :: numptr_amode(:)
    integer,  intent(in)    :: numptrcw_amode(:)
    real(r8), intent(in)    :: pi
    ! accum-coarse-exchange path flags (host-owned; passed in)
    integer,  intent(in)    :: modeptr_accum
    integer,  intent(in)    :: modeptr_coarse
    integer,  intent(in)    :: modeptr_stracoar
    integer,  intent(in)    :: igrow_shrink_renamexf(:)
    integer,  intent(in)    :: ixferable_all_renamexf(:)
    integer,  intent(in)    :: ixferable_a_renamexf(:,:)
    integer,  intent(in)    :: ixferable_c_renamexf(:,:)
    logical,  intent(in)    :: strat_only_renamexf(:)
    logical,  intent(in)    :: modal_accum_coarse_exch  ! select accum-coarse exchange path
    integer,  intent(in)    :: pver                 ! number of vertical levels
    real(r8), intent(in)    :: gravit               ! gravitational acceleration (m/s2)
    character(len=*), intent(out) :: errmsg         ! error message
    integer,  intent(out)   :: errflg               ! error flag

    errmsg = ''
    errflg = 0

    if (modal_accum_coarse_exch) then
       call modal_aero_rename_acc_crs_sub(        &
            ncol,                                   &
            loffset,           deltat,              &
            pdel,              troplev,             &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx,           &
            dqdt_rnpos,                             &
            ntot_amode,        npair_renamexf,      &
            modefrm_renamexf,  modetoo_renamexf,    &
            nspecfrm_renamexf,                      &
            lspecfrma_renamexf, lspecfrmc_renamexf, &
            lspectooa_renamexf, lspectooc_renamexf, &
            alnsg_amode,       voltonumblo_amode,   &
            voltonumbhi_amode, dgnum_amode,         &
            nspec_amode,       specmw_amode,        &
            specdens_amode,    lmassptr_amode,      &
            lmassptrcw_amode,  numptr_amode,        &
            numptrcw_amode,    pi,                  &
            modeptr_accum,     modeptr_coarse,      &
            modeptr_stracoar,                       &
            igrow_shrink_renamexf,                  &
            ixferable_all_renamexf,                 &
            ixferable_a_renamexf, ixferable_c_renamexf, &
            strat_only_renamexf,                    &
            pver,              gravit,              &
            errmsg,            errflg               )
    else
       ! no_acc path does not produce dqdt_rnpos; define the required output here.
       dqdt_rnpos(:,:,:) = 0.0_r8
       call modal_aero_rename_no_acc_crs_sub(             &
            ncol,                                   &
            loffset,           deltat,              &
            pdel,                                   &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx,           &
            ntot_amode,        npair_renamexf,      &
            modefrm_renamexf,  modetoo_renamexf,    &
            nspecfrm_renamexf,                      &
            lspecfrma_renamexf, lspecfrmc_renamexf, &
            lspectooa_renamexf, lspectooc_renamexf, &
            alnsg_amode,       voltonumblo_amode,   &
            voltonumbhi_amode, dgnum_amode,         &
            nspec_amode,       specmw_amode,        &
            specdens_amode,    lmassptr_amode,      &
            lmassptrcw_amode,  numptr_amode,        &
            numptrcw_amode,    pi,                  &
            pver,              gravit,              &
            errmsg,            errflg               )
    end if
  end subroutine modal_aero_rename_run

  !----------------------------------------------------------------------
  ! private methods

  subroutine modal_aero_rename_no_acc_crs_sub(                       &
                        ncol,                                   &
                        loffset,           deltat,              &
                        pdel,                                   &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx,              &
                        ntot_amode,        npair_renamexf,      &
                        modefrm_renamexf,  modetoo_renamexf,    &
                        nspecfrm_renamexf,                      &
                        lspecfrma_renamexf, lspecfrmc_renamexf, &
                        lspectooa_renamexf, lspectooc_renamexf, &
                        alnsg_amode,       voltonumblo_amode,   &
                        voltonumbhi_amode, dgnum_amode,         &
                        nspec_amode,       specmw_amode,        &
                        specdens_amode,    lmassptr_amode,      &
                        lmassptrcw_amode,  numptr_amode,        &
                        numptrcw_amode,    pi,                  &
                        pver,              gravit,              &
                        errmsg,            errflg               )
   use shr_spfn_mod, only: erfc => shr_spfn_erfc

   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)

   real(r8), intent(in)    :: pdel(:,:)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: q(:,:,:) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: qqcw(:,:,:) ! like q but for cloud-borne species

   real(r8), intent(inout) :: dqdt(:,:,:)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the
                              !     "fromwhere" continuous growth process
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(inout) :: dqqcwdt(:,:,:)
   real(r8), intent(in)    :: dqdt_other(:,:,:)
                              ! tendencies for "other" continuous growth process
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: dqqcwdt_other(:,:,:)
   logical,  intent(inout) :: dotendrn(:) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(:)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(:,:) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), intent(out)   :: qsrflx(:,:,:)
                              ! process-specific column tracer tendencies
   real(r8), intent(out)   :: qqcwsrflx(:,:,:)

   integer,  intent(in)    :: pver                 ! number of vertical levels
   real(r8), intent(in)    :: gravit               ! gravitational acceleration (m/s2)
   character(len=*), intent(out) :: errmsg         ! error message
   integer,  intent(out)   :: errflg               ! error flag
   ! shared mode metadata + resolved renaming-pair tables (host-owned; passed in)
   integer,  intent(in)    :: ntot_amode           ! number of aerosol modes
   integer,  intent(in)    :: npair_renamexf       ! number of renaming pairs
   integer,  intent(in)    :: modefrm_renamexf(:)  ! source mode index per pair
   integer,  intent(in)    :: modetoo_renamexf(:)  ! destination mode index per pair
   integer,  intent(in)    :: nspecfrm_renamexf(:) ! number of transferred species per pair
   integer,  intent(in)    :: lspecfrma_renamexf(:,:) ! interstitial source species (pcnst-space)
   integer,  intent(in)    :: lspecfrmc_renamexf(:,:) ! cloud-borne source species
   integer,  intent(in)    :: lspectooa_renamexf(:,:) ! interstitial destination species
   integer,  intent(in)    :: lspectooc_renamexf(:,:) ! cloud-borne destination species
   real(r8), intent(in)    :: alnsg_amode(:)       ! ln(geometric std dev) of each mode
   real(r8), intent(in)    :: voltonumblo_amode(:) ! volume-to-number ratio, low limit
   real(r8), intent(in)    :: voltonumbhi_amode(:) ! volume-to-number ratio, high limit
   real(r8), intent(in)    :: dgnum_amode(:)       ! nominal geometric mean diameter
   integer,  intent(in)    :: nspec_amode(:)       ! number of species in each mode
   real(r8), intent(in)    :: specmw_amode(:,:)    ! species molecular weight
   real(r8), intent(in)    :: specdens_amode(:,:)  ! species density
   integer,  intent(in)    :: lmassptr_amode(:,:)  ! interstitial mass pointer (pcnst-space)
   integer,  intent(in)    :: lmassptrcw_amode(:,:) ! cloud-borne mass pointer
   integer,  intent(in)    :: numptr_amode(:)      ! interstitial number pointer
   integer,  intent(in)    :: numptrcw_amode(:)    ! cloud-borne number pointer
   real(r8), intent(in)    :: pi                   ! pi
! !DESCRIPTION:
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1=-1
   integer :: i, icol_diag, ipair, iq, j, k, l, l1, la, lc, lunout
   integer :: lsfrma, lsfrmc, lstooa, lstooc
   integer :: mfrm, mtoo, n, n1, n2, ntot_msa_a
   integer :: idomode(ntot_amode)

   real (r8) :: deldryvol_a(ncol,pver,ntot_amode)
   real (r8) :: deldryvol_c(ncol,pver,ntot_amode)
   real (r8) :: deltatinv
   real (r8) :: dp_belowcut(maxpair_renamexf)
   real (r8) :: dp_cut(maxpair_renamexf)
   real (r8) :: dgn_aftr, dgn_xfer
   real (r8) :: dgn_t_new, dgn_t_old
   real (r8) :: dryvol_t_del, dryvol_t_new
   real (r8) :: dryvol_t_old, dryvol_t_oldbnd
   real (r8) :: dryvol_a(ncol,pver,ntot_amode)
   real (r8) :: dryvol_c(ncol,pver,ntot_amode)
   real (r8) :: dryvol_smallest(ntot_amode)
   real (r8) :: dum
   real (r8) :: dum3alnsg2(maxpair_renamexf)
   real (r8) :: dum_m2v, dum_m2vdt
   real (r8) :: factoraa(ntot_amode)
   real (r8) :: factoryy(ntot_amode)
   real (r8) :: frelax
   real (r8) :: lndp_cut(maxpair_renamexf)
   real (r8) :: lndgn_new, lndgn_old
   real (r8) :: lndgv_new, lndgv_old
   real (r8) :: num_t_old, num_t_oldbnd
   real (r8) :: onethird
   real (r8) :: pdel_fac
   real (r8) :: tailfr_volnew, tailfr_volold
   real (r8) :: tailfr_numnew, tailfr_numold
   real (r8) :: v2nhirlx(ntot_amode), v2nlorlx(ntot_amode)
   real (r8) :: xfercoef, xfertend
   real (r8) :: xferfrac_vol, xferfrac_num, xferfrac_max

   real (r8) :: yn_tail, yv_tail

! begin
  lunout = iulog
  errmsg = ''
  errflg = 0
  ! intent(out): fully define before any early return
  qsrflx(:,:,:) = 0.0_r8
  qqcwsrflx(:,:,:) = 0.0_r8

!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
! if (npair_renamexf .eq. -123456789) then
!     npair_renamexf = 0
!     call modal_aero_rename_init
! end if

!
!   check if any renaming pairs exist
!
  if (npair_renamexf .le. 0) return
!   if (ncol .ne. -123456789) return
! if (fromwhere .eq. 'aqchem') return

!
!   compute aerosol dry-volume for the "from mode" of each renaming pair
!   also compute dry-volume change during the continuous growth process
! using the incoming dqdt*deltat
!
  deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
  onethird = 1.0_r8/3.0_r8
  frelax = 27.0_r8
  xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

  do n = 1, ntot_amode
      idomode(n) = 0
  end do

  do ipair = 1, npair_renamexf
      if (ipair .gt. 1) goto 8100
      idomode(modefrm_renamexf(ipair)) = 1

      mfrm = modefrm_renamexf(ipair)
      mtoo = modetoo_renamexf(ipair)
      factoraa(mfrm) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mfrm)**2))
      factoraa(mtoo) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mtoo)**2))
      factoryy(mfrm) = sqrt( 0.5_r8 )/alnsg_amode(mfrm)
!   dryvol_smallest is a very small volume mixing ratio (m3-AP/kmol-air)
!   used for avoiding overflow.  it corresponds to dp = 1 nm
!   and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air
      dryvol_smallest(mfrm) = 1.0e-25_r8
      v2nlorlx(mfrm) = voltonumblo_amode(mfrm)*frelax
      v2nhirlx(mfrm) = voltonumbhi_amode(mfrm)/frelax

      dum3alnsg2(ipair) = 3.0_r8 * (alnsg_amode(mfrm)**2)
      dp_cut(ipair) = sqrt(   &
    dgnum_amode(mfrm)*exp(1.5_r8*(alnsg_amode(mfrm)**2)) *   &
    dgnum_amode(mtoo)*exp(1.5_r8*(alnsg_amode(mtoo)**2)) )
      lndp_cut(ipair) = log( dp_cut(ipair) )
      dp_belowcut(ipair) = 0.99_r8*dp_cut(ipair)
  end do

  do n = 1, ntot_amode
      if (idomode(n) .gt. 0) then
    dryvol_a(1:ncol,:,n) = 0.0_r8
    dryvol_c(1:ncol,:,n) = 0.0_r8
    deldryvol_a(1:ncol,:,n) = 0.0_r8
    deldryvol_c(1:ncol,:,n) = 0.0_r8
    do l1 = 1, nspec_amode(n)
!   dum_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!            [m3-AP/kmol-AP]= [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
                    dum_m2v = specmw_amode(l1,n) / specdens_amode(l1,n)
        dum_m2vdt = dum_m2v*deltat
        la = lmassptr_amode(l1,n)-loffset
        if (la > 0) then
        dryvol_a(1:ncol,:,n) = dryvol_a(1:ncol,:,n)    &
      + dum_m2v*max( 0.0_r8,   &
                          q(1:ncol,:,la)-deltat*dqdt_other(1:ncol,:,la) )
        deldryvol_a(1:ncol,:,n) = deldryvol_a(1:ncol,:,n)    &
      + (dqdt_other(1:ncol,:,la) + dqdt(1:ncol,:,la))*dum_m2vdt
        end if

        lc = lmassptrcw_amode(l1,n)-loffset
        if (lc > 0) then
        dryvol_c(1:ncol,:,n) = dryvol_c(1:ncol,:,n)    &
      + dum_m2v*max( 0.0_r8,   &
                          qqcw(1:ncol,:,lc)-deltat*dqqcwdt_other(1:ncol,:,lc) )
        deldryvol_c(1:ncol,:,n) = deldryvol_c(1:ncol,:,n)    &
      + (dqqcwdt_other(1:ncol,:,lc) +   &
               dqqcwdt(1:ncol,:,lc))*dum_m2vdt
        end if
    end do
      end if
  end do



!
!   loop over levels and columns to calc the renaming
!
mainloop1_k:  do k = 1, pver
mainloop1_i:  do i = 1, ncol

!   if dorename_atik is provided, then check if renaming needed at this i,k
  if (is_dorename_atik) then
      if (.not. dorename_atik(i,k)) cycle mainloop1_i
  end if
  pdel_fac = pdel(i,k)/gravit

!
!   loop over renameing pairs
!
mainloop1_ipair:  do ipair = 1, npair_renamexf

  mfrm = modefrm_renamexf(ipair)
  mtoo = modetoo_renamexf(ipair)

!   dryvol_t_old is the old total (a+c) dry-volume for the "from" mode
! in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
! (old/new = before/after the continuous growth)
  dryvol_t_old = dryvol_a(i,k,mfrm) + dryvol_c(i,k,mfrm)
  dryvol_t_del = deldryvol_a(i,k,mfrm) + deldryvol_c(i,k,mfrm)
  dryvol_t_new = dryvol_t_old + dryvol_t_del
  dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )

!   no renaming if dryvol_t_new ~ 0 or dryvol_t_del ~ 0
  if (dryvol_t_new .le. dryvol_smallest(mfrm)) cycle mainloop1_ipair
  if (dryvol_t_del .le. 1.0e-6_r8*dryvol_t_oldbnd) cycle mainloop1_ipair

!   num_t_old is total number in particles/kmol-air
  num_t_old = q(i,k,numptr_amode(mfrm)-loffset)
  num_t_old = num_t_old + qqcw(i,k,numptrcw_amode(mfrm)-loffset)
  num_t_old = max( 0.0_r8, num_t_old )
  dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
  num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
  num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   no renaming if dgnum < "base" dgnum,
  dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
  if (dgn_t_new .le. dgnum_amode(mfrm)) cycle mainloop1_ipair

!   compute new fraction of number and mass in the tail (dp > dp_cut)
  lndgn_new = log( dgn_t_new )
  lndgv_new = lndgn_new + dum3alnsg2(ipair)
  yn_tail = (lndp_cut(ipair) - lndgn_new)*factoryy(mfrm)
  yv_tail = (lndp_cut(ipair) - lndgv_new)*factoryy(mfrm)
  tailfr_numnew = 0.5_r8*erfc( yn_tail )
  tailfr_volnew = 0.5_r8*erfc( yv_tail )

!   compute old fraction of number and mass in the tail (dp > dp_cut)
  dgn_t_old =   &
    (dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
!   if dgn_t_new exceeds dp_cut, use the minimum of dgn_t_old and
!   dp_belowcut to guarantee some transfer
  if (dgn_t_new .ge. dp_cut(ipair)) then
      dgn_t_old = min( dgn_t_old, dp_belowcut(ipair) )
  end if
  lndgn_old = log( dgn_t_old )
  lndgv_old = lndgn_old + dum3alnsg2(ipair)
  yn_tail = (lndp_cut(ipair) - lndgn_old)*factoryy(mfrm)
  yv_tail = (lndp_cut(ipair) - lndgv_old)*factoryy(mfrm)
  tailfr_numold = 0.5_r8*erfc( yn_tail )
  tailfr_volold = 0.5_r8*erfc( yv_tail )

!   transfer fraction is difference between new and old tail-fractions
!   transfer fraction for number cannot exceed that of mass
  dum = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_old
  if (dum .le. 0.0_r8) cycle mainloop1_ipair

  xferfrac_vol = min( dum, dryvol_t_new )/dryvol_t_new
  xferfrac_vol = min( xferfrac_vol, xferfrac_max )
  xferfrac_num = tailfr_numnew - tailfr_numold
  xferfrac_num = max( 0.0_r8, min( xferfrac_num, xferfrac_vol ) )

!
!   compute tendencies for the renaming transfer
!
  j = jsrflx_rename
  do iq = 1, nspecfrm_renamexf(ipair)
      xfercoef = xferfrac_vol*deltatinv
      if (iq .eq. 1) xfercoef = xferfrac_num*deltatinv

      lsfrma = lspecfrma_renamexf(iq,ipair)-loffset
      lsfrmc = lspecfrmc_renamexf(iq,ipair)-loffset
      lstooa = lspectooa_renamexf(iq,ipair)-loffset
      lstooc = lspectooc_renamexf(iq,ipair)-loffset

      if (lsfrma .gt. 0) then
    xfertend = xfercoef*max( 0.0_r8,   &
          (q(i,k,lsfrma)+dqdt(i,k,lsfrma)*deltat) )

!   diagnostic output start ----------------------------------------
                if (ldiag1 > 0) then
                if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
                  if (lstooa .gt. 0) then
                    write(*,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend),   &
                        deltat*dqdt(i,k,lstooa), deltat*(dqdt(i,k,lstooa) + xfertend)
                  else
                    write(*,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend)
                  end if
                end if
                end if
!   diagnostic output end   ------------------------------------------


    dqdt(i,k,lsfrma) = dqdt(i,k,lsfrma) - xfertend
    qsrflx(i,lsfrma,j) = qsrflx(i,lsfrma,j) - xfertend*pdel_fac
    if (lstooa .gt. 0) then
        dqdt(i,k,lstooa) = dqdt(i,k,lstooa) + xfertend
        qsrflx(i,lstooa,j) = qsrflx(i,lstooa,j) + xfertend*pdel_fac
    end if
      end if

      if (lsfrmc .gt. 0) then
    xfertend = xfercoef*max( 0.0_r8,   &
          (qqcw(i,k,lsfrmc)+dqqcwdt(i,k,lsfrmc)*deltat) )
    dqqcwdt(i,k,lsfrmc) = dqqcwdt(i,k,lsfrmc) - xfertend
    qqcwsrflx(i,lsfrmc,j) = qqcwsrflx(i,lsfrmc,j) - xfertend*pdel_fac
    if (lstooc .gt. 0) then
        dqqcwdt(i,k,lstooc) = dqqcwdt(i,k,lstooc) + xfertend
        qqcwsrflx(i,lstooc,j) = qqcwsrflx(i,lstooc,j) + xfertend*pdel_fac
    end if
      end if

  end do   ! "iq = 1, nspecfrm_renamexf(ipair)"


  end do mainloop1_ipair


  end do mainloop1_i
  end do mainloop1_k

!
!   set dotend's
!
  dotendrn(:) = .false.
  dotendqqcwrn(:) = .false.
  do ipair = 1, npair_renamexf
  do iq = 1, nspecfrm_renamexf(ipair)
      lsfrma = lspecfrma_renamexf(iq,ipair) - loffset
      lsfrmc = lspecfrmc_renamexf(iq,ipair) - loffset
      lstooa = lspectooa_renamexf(iq,ipair) - loffset
      lstooc = lspectooc_renamexf(iq,ipair) - loffset
      if (lsfrma .gt. 0) then
    dotendrn(lsfrma) = .true.
    if (lstooa .gt. 0) dotendrn(lstooa) = .true.
      end if
      if (lsfrmc .gt. 0) then
    dotendqqcwrn(lsfrmc) = .true.
    if (lstooc .gt. 0) dotendqqcwrn(lstooc) = .true.
      end if
  end do
  end do


  return


!
!   error -- renaming currently just works for 1 pair
!
8100  write(lunout,9050) ipair
  errflg = 1
  errmsg = 'modal_aero_rename_no_acc_crs_sub error'
  return
9050  format( / '*** subr. modal_aero_rename_no_acc_crs_sub ***' /   &
            4x, 'aerosol renaming not implemented for ipair =', i5 )

!EOC
  end subroutine modal_aero_rename_no_acc_crs_sub



  subroutine modal_aero_rename_acc_crs_sub(                       &
                        ncol,                                   &
                        loffset,           deltat,              &
                        pdel,              troplev,             &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx,           &
                        dqdt_rnpos,              &
                        ntot_amode,        npair_renamexf,      &
                        modefrm_renamexf,  modetoo_renamexf,    &
                        nspecfrm_renamexf,                      &
                        lspecfrma_renamexf, lspecfrmc_renamexf, &
                        lspectooa_renamexf, lspectooc_renamexf, &
                        alnsg_amode,       voltonumblo_amode,   &
                        voltonumbhi_amode, dgnum_amode,         &
                        nspec_amode,       specmw_amode,        &
                        specdens_amode,    lmassptr_amode,      &
                        lmassptrcw_amode,  numptr_amode,        &
                        numptrcw_amode,    pi,                  &
                        modeptr_accum,     modeptr_coarse,      &
                        modeptr_stracoar,                       &
                        igrow_shrink_renamexf,                  &
                        ixferable_all_renamexf,                 &
                        ixferable_a_renamexf, ixferable_c_renamexf, &
                        strat_only_renamexf,                    &
                        pver,              gravit,              &
                        errmsg,            errflg               )

! !USES:

   use shr_spfn_mod, only: erfc => shr_spfn_erfc

! !PARAMETERS:
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)
   integer,  intent(in)    :: troplev(:)

   real(r8), intent(in)    :: pdel(:,:)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: q(:,:,:) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: qqcw(:,:,:) ! like q but for cloud-borne species

   real(r8), intent(inout) :: dqdt(:,:,:)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the
                              !     "fromwhere" continuous growth process
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(inout) :: dqqcwdt(:,:,:)
   real(r8), intent(in)    :: dqdt_other(:,:,:)
                              ! tendencies for "other" continuous growth process
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: dqqcwdt_other(:,:,:)
   logical,  intent(inout) :: dotendrn(:) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(:)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(:,:) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), intent(out)   :: qsrflx(:,:,:)
                              ! process-specific column tracer tendencies
   real(r8), intent(out)   :: qqcwsrflx(:,:,:)
   real(r8), intent(out)   :: dqdt_rnpos(:,:,:)
                              ! the positive (production) part of the renaming tendency

   integer,  intent(in)    :: pver                 ! number of vertical levels
   real(r8), intent(in)    :: gravit               ! gravitational acceleration (m/s2)
   character(len=*), intent(out) :: errmsg         ! error message
   integer,  intent(out)   :: errflg               ! error flag
   ! shared mode metadata + resolved renaming-pair tables (host-owned; passed in)
   integer,  intent(in)    :: ntot_amode           ! number of aerosol modes
   integer,  intent(in)    :: npair_renamexf       ! number of renaming pairs
   integer,  intent(in)    :: modefrm_renamexf(:)  ! source mode index per pair
   integer,  intent(in)    :: modetoo_renamexf(:)  ! destination mode index per pair
   integer,  intent(in)    :: nspecfrm_renamexf(:) ! number of transferred species per pair
   integer,  intent(in)    :: lspecfrma_renamexf(:,:) ! interstitial source species (pcnst-space)
   integer,  intent(in)    :: lspecfrmc_renamexf(:,:) ! cloud-borne source species
   integer,  intent(in)    :: lspectooa_renamexf(:,:) ! interstitial destination species
   integer,  intent(in)    :: lspectooc_renamexf(:,:) ! cloud-borne destination species
   real(r8), intent(in)    :: alnsg_amode(:)       ! ln(geometric std dev) of each mode
   real(r8), intent(in)    :: voltonumblo_amode(:) ! volume-to-number ratio, low limit
   real(r8), intent(in)    :: voltonumbhi_amode(:) ! volume-to-number ratio, high limit
   real(r8), intent(in)    :: dgnum_amode(:)       ! nominal geometric mean diameter
   integer,  intent(in)    :: nspec_amode(:)       ! number of species in each mode
   real(r8), intent(in)    :: specmw_amode(:,:)    ! species molecular weight
   real(r8), intent(in)    :: specdens_amode(:,:)  ! species density
   integer,  intent(in)    :: lmassptr_amode(:,:)  ! interstitial mass pointer (pcnst-space)
   integer,  intent(in)    :: lmassptrcw_amode(:,:) ! cloud-borne mass pointer
   integer,  intent(in)    :: numptr_amode(:)      ! interstitial number pointer
   integer,  intent(in)    :: numptrcw_amode(:)    ! cloud-borne number pointer
   real(r8), intent(in)    :: pi                   ! pi
   ! accum-coarse-exchange path flags (host-owned; passed in)
   integer,  intent(in)    :: modeptr_accum        ! accumulation mode index
   integer,  intent(in)    :: modeptr_coarse       ! coarse mode index
   integer,  intent(in)    :: modeptr_stracoar     ! stratospheric coarse mode index
   integer,  intent(in)    :: igrow_shrink_renamexf(:)  ! +1 growing / -1 shrinking per pair
   integer,  intent(in)    :: ixferable_all_renamexf(:) ! all-species-transferable flag per pair
   integer,  intent(in)    :: ixferable_a_renamexf(:,:) ! per-species interstitial transferable flag
   integer,  intent(in)    :: ixferable_c_renamexf(:,:) ! per-species cloud-borne transferable flag
   logical,  intent(in)    :: strat_only_renamexf(:)    ! restrict renaming to the stratosphere
! !DESCRIPTION:
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1 = -1
   integer :: i, icol_diag, ipair, iq
   integer :: j, k
   integer :: l, l1, la, lc, lunout
   integer :: lsfrma, lsfrmc, lstooa, lstooc
   integer :: mfrm, mtoo, n, n1, n2, ntot_msa_a

   logical :: l_dqdt_rnpos
   logical :: flagaa_shrink, flagbb_shrink

   real (r8) :: deldryvol_a(ncol,pver)
   real (r8) :: deldryvol_c(ncol,pver)
   real (r8) :: deltatinv
   real (r8) :: dgn_aftr, dgn_xfer
   real (r8) :: dgn_t_new, dgn_t_old, dgn_t_oldb
   real (r8) :: dryvol_t_del, dryvol_t_new, dryvol_t_new_xfab
   real (r8) :: dryvol_t_old, dryvol_t_oldb, dryvol_t_oldbnd
   real (r8) :: dryvol_a(ncol,pver)
   real (r8) :: dryvol_c(ncol,pver)
   real (r8) :: dryvol_a_xfab(ncol,pver)
   real (r8) :: dryvol_c_xfab(ncol,pver)
   real (r8) :: dryvol_xferamt
   real (r8) :: lndgn_new, lndgn_old
   real (r8) :: lndgv_new, lndgv_old
   real (r8) :: num_t_old, num_t_oldbnd
   real (r8) :: onethird
   real (r8) :: pdel_fac
   real (r8) :: tailfr_volnew, tailfr_volold
   real (r8) :: tailfr_numnew, tailfr_numold
   real (r8) :: tmpa, tmpf
   real (r8) :: tmp_m2v, tmp_m2vdt
   real (r8) :: xfercoef, xfertend
   real (r8) :: xferfrac_vol, xferfrac_num, xferfrac_max

   real (r8) :: yn_tail, yv_tail

! begin
  lunout = iulog
  errmsg = ''
  errflg = 0
  ! intent(out): fully define before any early return
  qsrflx(:,:,:) = 0.0_r8
  qqcwsrflx(:,:,:) = 0.0_r8

!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
! if (npair_renamexf .eq. -123456789) then
!     npair_renamexf = 0
!     call modal_aero_rename_init
! end if

!
!   check if any renaming pairs exist
!
  if (npair_renamexf .le. 0) return
!   if (ncol .ne. -123456789) return
! if (fromwhere .eq. 'aqchem') return


  deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
  onethird = 1.0_r8/3.0_r8
  xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

  ! dqdt_rnpos is now a required output; always produced.
  l_dqdt_rnpos = .true.
  dqdt_rnpos(:,:,:) = 0.0_r8



!
!   loop over renaming pairs
!
mainloop1_ipair:  do ipair = 1, npair_renamexf

  mfrm = modefrm_renamexf(ipair)
  mtoo = modetoo_renamexf(ipair)

        flagaa_shrink = &
            ((mfrm==modeptr_coarse) .and.  (mtoo==modeptr_accum)) .or. &
            ((mfrm==modeptr_stracoar) .and. (mtoo==modeptr_accum))

!
!   compute aerosol dry-volume for the "from mode" of each renaming pair
!   also compute dry-volume change during the continuous growth process
! using the incoming dqdt*deltat
!
  dryvol_a(:,:) = 0.0_r8
  dryvol_c(:,:) = 0.0_r8
  deldryvol_a(:,:) = 0.0_r8
  deldryvol_c(:,:) = 0.0_r8
  if (ixferable_all_renamexf(ipair) <= 0) then
      dryvol_a_xfab(:,:) = 0.0_r8
      dryvol_c_xfab(:,:) = 0.0_r8
  end if

  n = mfrm
  do l1 = 1, nspec_amode(n)
!   tmp_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!            [m3-AP/kmol-AP]= [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
            tmp_m2v = specmw_amode(l1,n) / specdens_amode(l1,n)
      tmp_m2vdt = tmp_m2v*deltat
      la = lmassptr_amode(l1,n)-loffset
      if (la > 0) then
    dryvol_a(1:ncol,:) = dryvol_a(1:ncol,:)    &
        + tmp_m2v*max( 0.0_r8,   &
          q(1:ncol,:,la)-deltat*dqdt_other(1:ncol,:,la) )
    deldryvol_a(1:ncol,:) = deldryvol_a(1:ncol,:)    &
        + (dqdt_other(1:ncol,:,la) + dqdt(1:ncol,:,la))*tmp_m2vdt
    if ( (ixferable_all_renamexf(ipair) <= 0) .and. &
         (ixferable_a_renamexf(l1,ipair) > 0) ) then
        dryvol_a_xfab(1:ncol,:) = dryvol_a_xfab(1:ncol,:)    &
      + tmp_m2v*max( 0.0_r8,   &
      q(1:ncol,:,la)+deltat*dqdt(1:ncol,:,la) )
    end if
      end if

      lc = lmassptrcw_amode(l1,n)-loffset
      if (lc > 0) then
    dryvol_c(1:ncol,:) = dryvol_c(1:ncol,:)    &
        + tmp_m2v*max( 0.0_r8,   &
          qqcw(1:ncol,:,lc)-deltat*dqqcwdt_other(1:ncol,:,lc) )
    deldryvol_c(1:ncol,:) = deldryvol_c(1:ncol,:)    &
        + (dqqcwdt_other(1:ncol,:,lc) +   &
                 dqqcwdt(1:ncol,:,lc))*tmp_m2vdt
    if ( (ixferable_all_renamexf(ipair) <= 0) .and. &
         (ixferable_c_renamexf(l1,ipair) > 0) ) then
        dryvol_c_xfab(1:ncol,:) = dryvol_c_xfab(1:ncol,:)    &
      + tmp_m2v*max( 0.0_r8,   &
        qqcw(1:ncol,:,lc)+deltat*dqqcwdt(1:ncol,:,lc) )
    end if
      end if
  end do

!
!
!   loop over levels and columns to calc the renaming
!
!
mainloop1_k:  do k = 1, pver
mainloop1_i:  do i = 1, ncol

!   if dorename_atik is provided, then check if renaming needed at this i,k
  if (is_dorename_atik) then
      if (.not. dorename_atik(i,k)) cycle mainloop1_i
  end if

!   if strat_only_renamexf is true, then cycle when at or below the tropopause level
        if ( strat_only_renamexf(ipair) ) then
            if ( k >= troplev(i) ) cycle mainloop1_i
        end if


!   dryvol_t_old is the old total (a+c) dry-volume for the "from" mode
! in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
! (old/new = before/after the continuous growth)
  dryvol_t_old = dryvol_a(i,k) + dryvol_c(i,k)
  dryvol_t_del = deldryvol_a(i,k) + deldryvol_c(i,k)
  dryvol_t_new = dryvol_t_old + dryvol_t_del
  dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )

grow_shrink_conditional1: &
  if (igrow_shrink_renamexf(ipair) > 0) then
!   do renaming for growing particles

!   no renaming if dryvol_t_new ~ 0
  if (dryvol_t_new .le. dryvol_smallest(mfrm)) cycle mainloop1_i
!   no renaming if delta_dryvol is very small or negative
  if ( (method_optbb_renamexf /= 2) .and. &
       (dryvol_t_del .le. 1.0e-6_r8*dryvol_t_oldbnd) ) cycle mainloop1_i

!   num_t_old is total number in particles/kmol-air
  num_t_old = q(i,k,numptr_amode(mfrm)-loffset)
  num_t_old = num_t_old + qqcw(i,k,numptrcw_amode(mfrm)-loffset)
  num_t_old = max( 0.0_r8, num_t_old )
  dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
  num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
  num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   compute new dgnum
  dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
!   no renaming if dgn_t_new < threshold value
  if (dgn_t_new .le. dp_xfernone_threshaa(ipair)) cycle mainloop1_i

!   compute old dgnum and possibly a smaller value to get more renaming transfer
  dgn_t_old =   &
    (dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
  dgn_t_oldb = dgn_t_old
  dryvol_t_oldb = dryvol_t_old
  if ( method_optbb_renamexf == 2) then
      if (dgn_t_old .ge. dp_cut(ipair)) then
    ! this revised volume corresponds to dgn_t_old == dp_belowcut, and same number conc
    dryvol_t_oldb = dryvol_t_old * (dp_belowcut(ipair)/dgn_t_old)**3
    dgn_t_oldb = dp_belowcut(ipair)
      end if
      if (dgn_t_new .lt. dp_xferall_thresh(ipair)) then
    !   no renaming if delta_dryvol is very small or negative
    if ((dryvol_t_new-dryvol_t_oldb) .le. 1.0e-6_r8*dryvol_t_oldbnd) cycle mainloop1_i
      end if

  else if (dgn_t_new .ge. dp_cut(ipair)) then
!   if dgn_t_new exceeds dp_cut, use the minimum of dgn_t_oldb and
!   dp_belowcut to guarantee some transfer
      dgn_t_oldb = min( dgn_t_oldb, dp_belowcut(ipair) )
  end if

!   compute new fraction of number and mass in the tail (dp > dp_cut)
  lndgn_new = log( dgn_t_new )
  lndgv_new = lndgn_new + factor_3alnsg2(ipair)
  yn_tail = (lndp_cut(ipair) - lndgn_new)*factoryy(mfrm)
  yv_tail = (lndp_cut(ipair) - lndgv_new)*factoryy(mfrm)
  tailfr_numnew = 0.5_r8*erfc( yn_tail )
  tailfr_volnew = 0.5_r8*erfc( yv_tail )

!   compute old fraction of number and mass in the tail (dp > dp_cut)
  lndgn_old = log( dgn_t_oldb )
  lndgv_old = lndgn_old + factor_3alnsg2(ipair)
  yn_tail = (lndp_cut(ipair) - lndgn_old)*factoryy(mfrm)
  yv_tail = (lndp_cut(ipair) - lndgv_old)*factoryy(mfrm)
  tailfr_numold = 0.5_r8*erfc( yn_tail )
  tailfr_volold = 0.5_r8*erfc( yv_tail )

!   transfer fraction is difference between new and old tail-fractions
!   transfer fraction for number cannot exceed that of mass
  if ( (method_optbb_renamexf == 2) .and. &
       (dgn_t_new .ge. dp_xferall_thresh(ipair)) ) then
      dryvol_xferamt = dryvol_t_new
  else
      dryvol_xferamt = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_oldb
  end if
  if (dryvol_xferamt .le. 0.0_r8) cycle mainloop1_i

  xferfrac_vol = max( 0.0_r8, (dryvol_xferamt/dryvol_t_new) )
  if ( method_optbb_renamexf == 2 .and. &
       (xferfrac_vol >= xferfrac_max) ) then
      ! transfer entire contents of mode
      xferfrac_vol = 1.0_r8
      xferfrac_num = 1.0_r8
  else
      xferfrac_vol = min( xferfrac_vol, xferfrac_max )
      xferfrac_num = tailfr_numnew - tailfr_numold
      xferfrac_num = max( 0.0_r8, min( xferfrac_num, xferfrac_vol ) )
  end if

  if (ixferable_all_renamexf(ipair) <= 0) then
      ! not all species are xferable
      dryvol_t_new_xfab = max( 0.0_r8, (dryvol_a_xfab(i,k) + dryvol_c_xfab(i,k)) )
      dryvol_xferamt = xferfrac_vol*dryvol_t_new
      if (dryvol_t_new_xfab >= 0.999999_r8*dryvol_xferamt) then
    ! xferable dryvol can supply the needed dryvol_xferamt
    ! but xferfrac_vol must be increased
    xferfrac_vol = min( 1.0_r8, (dryvol_xferamt/dryvol_t_new_xfab) )
      else if (dryvol_t_new_xfab >= 1.0e-7_r8*dryvol_xferamt) then
    ! xferable dryvol cannot supply the needed dryvol_xferamt
    ! so transfer all of it, and reduce the number transfer
    xferfrac_vol = 1.0_r8
    xferfrac_num = xferfrac_num*(dryvol_t_new_xfab/dryvol_xferamt)
      else
    ! xferable dryvol << needed dryvol_xferamt
    cycle mainloop1_i
      end if
  end if

  else grow_shrink_conditional1
!   do renaming for shrinking particles

!   no renaming if (dryvol_t_old ~ 0)
  if (dryvol_t_old .le. dryvol_smallest(mfrm)) cycle mainloop1_i

!   when (delta_dryvol is very small or positive),
!      which means particles are not evaporating,
!      only do renaming if [(flagaa_shrink true) and (in stratosphere)]],
!   and set flagbb_shrink true to identify this special case
  if (dryvol_t_del .ge. -1.0e-6_r8*dryvol_t_oldbnd) then
      if ( ( flagaa_shrink ) .and. ( k < troplev(i) ) ) then
    flagbb_shrink = .true.
      else
    cycle mainloop1_i
      end if
  else
      flagbb_shrink = .false.
  end if

!   num_t_old is total number in particles/kmol-air
  num_t_old = q(i,k,numptr_amode(mfrm)-loffset)
  num_t_old = num_t_old + qqcw(i,k,numptrcw_amode(mfrm)-loffset)
  num_t_old = max( 0.0_r8, num_t_old )
  dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
  num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
  num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   compute new dgnum
  dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
!   no renaming if (dgn_t_new > xfernone threshold value)
  if (dgn_t_new .ge. dp_xfernone_threshaa(ipair)) cycle mainloop1_i
!   if (flagbb_shrink true), renaming only when (dgn_t_new <= dp_cut value)
  if ( flagbb_shrink ) then
      if (dgn_t_new .gt. dp_cut(ipair)) cycle mainloop1_i
  end if

  if ( dgn_t_new .le. dp_xferall_thresh(ipair) ) then
!   special case of (dgn_t_new <= xferall threshold value)
      tailfr_numnew = 1.0_r8
      tailfr_volnew = 1.0_r8
  else
!   compute new fraction of number and mass in the tail (dp < dp_cut)
      lndgn_new = log( dgn_t_new )
      lndgv_new = lndgn_new + factor_3alnsg2(ipair)
      yn_tail = (lndp_cut(ipair) - lndgn_new)*factoryy(mfrm)
      yv_tail = (lndp_cut(ipair) - lndgv_new)*factoryy(mfrm)
      tailfr_numnew = 1.0_r8 - 0.5_r8*erfc( yn_tail )
      tailfr_volnew = 1.0_r8 - 0.5_r8*erfc( yv_tail )
  end if

!   compute old dgnum
  dgn_t_old =   &
    (dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
  dgn_t_oldb = dgn_t_old
  dryvol_t_oldb = dryvol_t_old

!   no need to compute old fraction of number and mass in the tail
  tailfr_numold = 0.0_r8
  tailfr_volold = 0.0_r8

!   transfer fraction is new tail-fraction
  xferfrac_vol = tailfr_volnew
  if (xferfrac_vol .le. 0.0_r8) cycle mainloop1_i
  xferfrac_num = tailfr_numnew

  if (xferfrac_vol >= xferfrac_max) then
      ! transfer entire contents of mode
      xferfrac_vol = 1.0_r8
      xferfrac_num = 1.0_r8
  else
      xferfrac_vol = min( xferfrac_vol, xferfrac_max )
!   transfer fraction for number cannot be less than that of volume
      xferfrac_num = max( xferfrac_num, xferfrac_vol )
      xferfrac_num = min( xferfrac_max, xferfrac_num )
  end if

  if (ixferable_all_renamexf(ipair) <= 0) then
      ! not all species are xferable
      dryvol_t_new_xfab = max( 0.0_r8, (dryvol_a_xfab(i,k) + dryvol_c_xfab(i,k)) )
      dryvol_xferamt = xferfrac_vol*dryvol_t_new
      if (dryvol_t_new_xfab >= 0.999999_r8*dryvol_xferamt) then
    ! xferable dryvol can supply the needed dryvol_xferamt
    ! but xferfrac_vol must be increased
    xferfrac_vol = min( 1.0_r8, (dryvol_xferamt/dryvol_t_new_xfab) )
      else if (dryvol_t_new_xfab >= 1.0e-7_r8*dryvol_xferamt) then
    ! xferable dryvol cannot supply the needed dryvol_xferamt
    ! so transfer all of it, and reduce the number transfer
    xferfrac_vol = 1.0_r8
    xferfrac_num = xferfrac_num*(dryvol_t_new_xfab/dryvol_xferamt)
      else
    ! xferable dryvol << needed dryvol_xferamt
    cycle mainloop1_i
      end if
  end if

  endif grow_shrink_conditional1

!
!   compute tendencies for the renaming transfer
!
  pdel_fac = pdel(i,k)/gravit
  j = jsrflx_rename
  do iq = 1, nspecfrm_renamexf(ipair)
      xfercoef = xferfrac_vol*deltatinv
      if (iq .eq. 1) xfercoef = xferfrac_num*deltatinv

      lsfrma = lspecfrma_renamexf(iq,ipair)-loffset
      lsfrmc = lspecfrmc_renamexf(iq,ipair)-loffset
      lstooa = lspectooa_renamexf(iq,ipair)-loffset
      lstooc = lspectooc_renamexf(iq,ipair)-loffset

      if (lsfrma .gt. 0) then
    xfertend = xfercoef*max( 0.0_r8,   &
          (q(i,k,lsfrma)+dqdt(i,k,lsfrma)*deltat) )

!   diagnostic output start ----------------------------------------
                if (ldiag1 > 0) then
                if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
                  if (lstooa .gt. 0) then
                    write(iulog,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend),   &
                        deltat*dqdt(i,k,lstooa), deltat*(dqdt(i,k,lstooa) + xfertend)
                  else
                    write(iulog,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend)
                  end if
                end if
                end if
!   diagnostic output end   ------------------------------------------


    dqdt(i,k,lsfrma) = dqdt(i,k,lsfrma) - xfertend
    qsrflx(i,lsfrma,j) = qsrflx(i,lsfrma,j) - xfertend*pdel_fac
    if (lstooa .gt. 0) then
        dqdt(i,k,lstooa) = dqdt(i,k,lstooa) + xfertend
        qsrflx(i,lstooa,j) = qsrflx(i,lstooa,j) + xfertend*pdel_fac
        if ( l_dqdt_rnpos ) &
      dqdt_rnpos(i,k,lstooa) = dqdt_rnpos(i,k,lstooa) + xfertend
    end if
      end if

      if (lsfrmc .gt. 0) then
    xfertend = xfercoef*max( 0.0_r8,   &
          (qqcw(i,k,lsfrmc)+dqqcwdt(i,k,lsfrmc)*deltat) )
    dqqcwdt(i,k,lsfrmc) = dqqcwdt(i,k,lsfrmc) - xfertend
    qqcwsrflx(i,lsfrmc,j) = qqcwsrflx(i,lsfrmc,j) - xfertend*pdel_fac
    if (lstooc .gt. 0) then
        dqqcwdt(i,k,lstooc) = dqqcwdt(i,k,lstooc) + xfertend
        qqcwsrflx(i,lstooc,j) = qqcwsrflx(i,lstooc,j) + xfertend*pdel_fac
    end if
      end if

  end do   ! "iq = 1, nspecfrm_renamexf(ipair)"


  end do mainloop1_i
  end do mainloop1_k


  end do mainloop1_ipair

!
!   set dotend's
!
  dotendrn(:) = .false.
  dotendqqcwrn(:) = .false.
  do ipair = 1, npair_renamexf
  do iq = 1, nspecfrm_renamexf(ipair)
      lsfrma = lspecfrma_renamexf(iq,ipair) - loffset
      lsfrmc = lspecfrmc_renamexf(iq,ipair) - loffset
      lstooa = lspectooa_renamexf(iq,ipair) - loffset
      lstooc = lspectooc_renamexf(iq,ipair) - loffset
      if (lsfrma .gt. 0) then
    dotendrn(lsfrma) = .true.
    if (lstooa .gt. 0) dotendrn(lstooa) = .true.
      end if
      if (lsfrmc .gt. 0) then
    dotendqqcwrn(lsfrmc) = .true.
    if (lstooc .gt. 0) dotendqqcwrn(lstooc) = .true.
      end if
  end do
  end do


  return


!
!   error -- renaming currently just works for 1 pair
!
8100  write(lunout,9050) ipair
  errflg = 1
  errmsg = 'modal_aero_rename_acc_crs_sub error'
  return
9050  format( / '*** subr. modal_aero_rename_acc_crs_sub ***' /   &
            4x, 'aerosol renaming not implemented for ipair =', i5 )

  end subroutine modal_aero_rename_acc_crs_sub
end module modal_aero_rename
