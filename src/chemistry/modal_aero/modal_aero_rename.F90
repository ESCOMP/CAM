! modal_aero_rename.F90
!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_rename --- modal aerosol mode merging (renaming)
!
! !INTERFACE:
  module modal_aero_rename

! !USES:
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use cam_abortutils,  only: endrun
  use cam_logfile,     only: iulog
  use mo_constants,    only: pi
  use chem_mods,       only: gas_pcnst
  use ppgrid,          only: pcols, pver
  use constituents,    only: pcnst, cnst_name
  use spmd_utils,      only: masterproc
  use modal_aero_data, only: maxspec_renamexf=>nspec_max, ntot_amode
  use modal_aero_data, only: alnsg_amode, voltonumblo_amode, voltonumbhi_amode, dgnum_amode, nspec_amode
  use modal_aero_data, only: specmw_amode, specdens_amode, lmassptr_amode, lmassptrcw_amode
  use modal_aero_data, only: numptr_amode, numptrcw_amode, modeptr_coarse, modeptr_accum
  use modal_aero_data, only: specmw_amode, specdens_amode, lmassptr_amode, lmassptrcw_amode, numptr_amode, numptrcw_amode
  use modal_aero_data, only: dgnumhi_amode, dgnumlo_amode, cnst_name_cw, modeptr_aitken

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public modal_aero_rename_sub, modal_aero_rename_init

! !PUBLIC DATA MEMBERS:
  integer, parameter :: pcnstxx = gas_pcnst

! *** select one of the 3 following options
! *** for maxpair_renamexf = 2 or 3, use mode definition files with
!     dgnumhi_amode(modeptr_accum)  = 1.1e-6 m
!     dgnumlo_amode(modeptr_coarse) = 0.9e-6 m

! integer, parameter, public :: maxpair_renamexf = 1
! integer, parameter, public :: ipair_select_renamexf(maxpair_renamexf) = (/ 2001 /)

! integer, parameter, public :: maxpair_renamexf = 2
! integer, parameter, public :: ipair_select_renamexf(maxpair_renamexf) = (/ 2001, 1003 /)

  integer, parameter, public :: maxpair_renamexf = 3
  integer, parameter, public :: ipair_select_renamexf(maxpair_renamexf) = (/ 2001, 1003, 3001 /)
! ipair_select_renamexf defines the mode_from and mode_too for each renaming pair
! 2001 = aitken --> accum
! 1003 = accum  --> coarse
! 3001 = coarse --> accum

  integer, parameter, public :: method_optbb_renamexf = 2

  integer, public :: npair_renamexf = -123456789
  integer, protected, public :: modefrm_renamexf(maxpair_renamexf)
  integer, protected, public :: modetoo_renamexf(maxpair_renamexf)
  integer, protected, public :: nspecfrm_renamexf(maxpair_renamexf)

  integer, allocatable, protected, public :: lspecfrma_renamexf(:,:)
  integer, allocatable, protected, public :: lspecfrmc_renamexf(:,:)
  integer, allocatable, protected, public :: lspectooa_renamexf(:,:)
  integer, allocatable, protected, public :: lspectooc_renamexf(:,:)

  integer, protected, public :: igrow_shrink_renamexf(maxpair_renamexf)
  integer, protected, public :: ixferable_all_renamexf(maxpair_renamexf)
  integer, protected, public :: ixferable_all_needed_renamexf(maxpair_renamexf)
  integer, allocatable, protected, public :: ixferable_a_renamexf(:,:)
  integer, allocatable, protected, public :: ixferable_c_renamexf(:,:)

  logical, public :: strat_only_renamexf(maxpair_renamexf)
! strat_only_renamexf - when true for a particular renaming pair, renaming is only
!                       done in stratosphere (when k < troplev(icol) )

! !PRIVATE DATA MEMBERS:
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

  logical :: modal_accum_coarse_exch = .false.

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

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_init(modal_accum_coarse_exch_in)
    logical, optional, intent(in) :: modal_accum_coarse_exch_in
    
    allocate( lspecfrma_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( lspecfrmc_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( lspectooa_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( lspectooc_renamexf(maxspec_renamexf,maxpair_renamexf) )

    allocate( ixferable_a_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( ixferable_c_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( ido_mode_calcaa(ntot_amode) )

    allocate( dryvol_smallest(ntot_amode) )
    allocate( factoraa(ntot_amode) )
    allocate( factoryy(ntot_amode) )

    allocate( v2nhirlx(ntot_amode), v2nlorlx(ntot_amode) )

    if (present(modal_accum_coarse_exch_in)) then
       modal_accum_coarse_exch = modal_accum_coarse_exch_in
    endif

    if (modal_accum_coarse_exch) then
       call modal_aero_rename_acc_crs_init()
    else
       call modal_aero_rename_no_acc_crs_init()
    endif

  end subroutine modal_aero_rename_init

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_sub(                       &
       fromwhere,         lchnk,               &
       ncol,              nstep,               &
       loffset,           deltat,              &
       pdel,              troplev,             &
       dotendrn,          q,                   &
       dqdt,              dqdt_other,          &
       dotendqqcwrn,      qqcw,                &
       dqqcwdt,           dqqcwdt_other,       &
       is_dorename_atik,  dorename_atik,       &
       jsrflx_rename,     nsrflx,              &
       qsrflx,            qqcwsrflx,           &
       dqdt_rnpos                              )


    ! !PARAMETERS:
    character(len=*), intent(in) :: fromwhere    ! identifies which module
    ! is making the call
    integer,  intent(in)    :: lchnk                ! chunk identifier
    integer,  intent(in)    :: ncol                 ! number of atmospheric column
    integer,  intent(in)    :: nstep                ! model time-step number
    integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
    real(r8), intent(in)    :: deltat               ! time step (s)
    integer,  intent(in)    :: troplev(pcols)

    real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
    real(r8), intent(in)    :: q(ncol,pver,pcnstxx) ! tracer mixing ratio array
    ! *** MUST BE mol/mol-air or #/mol-air
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(in)    :: qqcw(ncol,pver,pcnstxx) ! like q but for cloud-borne species

    real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)  ! TMR tendency array;
    ! incoming dqdt = tendencies for the 
    !     "fromwhere" continuous growth process 
    ! the renaming tendencies are added on
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx)  
    ! tendencies for "other" continuous growth process 
    ! currently in cam3
    !     dqdt is from gas (h2so4, nh3) condensation
    !     dqdt_other is from aqchem and soa
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx)  
    logical,  intent(inout) :: dotendrn(pcnstxx) ! identifies the species for which
    !     renaming dqdt is computed
    logical,  intent(inout) :: dotendqqcwrn(pcnstxx)

    logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
    logical,  intent(in)    :: dorename_atik(ncol,pver) ! true if renaming should
    ! be done at i,k
    integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
    integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

    real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    ! process-specific column tracer tendencies 
    real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    real(r8), optional, intent(out) &
         :: dqdt_rnpos(ncol,pver,pcnstxx)
    ! the positive (production) part of the renaming tendency

    if (modal_accum_coarse_exch) then
       call modal_aero_rename_acc_crs_sub(        &
            fromwhere,         lchnk,               &
            ncol,              nstep,               &
            loffset,           deltat,              &
            pdel,              troplev,             &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx,           &
            dqdt_rnpos                              )
    else
       call modal_aero_rename_no_acc_crs_sub(             &
            fromwhere,         lchnk,               &
            ncol,              nstep,               &
            loffset,           deltat,              &
            pdel,                                   &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx            )
    endif
  end subroutine modal_aero_rename_sub

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! private methods
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_rename_no_acc_crs_sub --- ...
!
! !INTERFACE:
	subroutine modal_aero_rename_no_acc_crs_sub(                       &
                        fromwhere,         lchnk,               &
                        ncol,              nstep,               &
                        loffset,           deltat,              &
                        pdel,                                   &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx            )

! !USES:
   use physconst, only: gravit, mwdry
   use units, only: getunit
   use shr_spfn_mod, only: erfc => shr_spfn_erfc

   implicit none


! !PARAMETERS:
   character(len=*), intent(in) :: fromwhere    ! identifies which module
                                                ! is making the call
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)

   real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: q(ncol,pver,pcnstxx) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: qqcw(ncol,pver,pcnstxx) ! like q but for cloud-borne species

   real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the 
                              !     "fromwhere" continuous growth process 
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
   real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx)  
                              ! tendencies for "other" continuous growth process 
                              ! currently in cam3
                              !     dqdt is from gas (h2so4, nh3) condensation
                              !     dqdt_other is from aqchem and soa
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx)  
   logical,  intent(inout) :: dotendrn(pcnstxx) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(pcnstxx)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(ncol,pver) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
                              ! process-specific column tracer tendencies 
   real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)

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
   integer, save :: lun = -1  ! logical unit for diagnostics (6, or other
                              ! if a special diagnostics file is opened)


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

!   get logical unit (for output to dumpconv, deactivate the "lun = 6")
 	lun = iulog
	if (lun < 1) then
	   lun = getunit()
 	   open( unit=lun, file='dump.rename',   &
 			status='unknown', form='formatted' )
	end if


!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
!	if (npair_renamexf .eq. -123456789) then
!	    npair_renamexf = 0
!	    call modal_aero_rename_init
!	end if

!
!   check if any renaming pairs exist
!
	if (npair_renamexf .le. 0) return
! 	if (ncol .ne. -123456789) return
!	if (fromwhere .eq. 'aqchem') return

!
!   compute aerosol dry-volume for the "from mode" of each renaming pair
!   also compute dry-volume change during the continuous growth process
!	using the incoming dqdt*deltat
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
!	in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
!	(old/new = before/after the continuous growth)
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

!   diagnostic output start ----------------------------------------
!!$ 	if (ldiag1 > 0) then
!!$ 	icol_diag = -1
!!$ 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) icol_diag = i
!!$ 	if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
!!$ !	write(lun,97010) fromwhere, nstep, lchnk, i, k, ipair
!!$ 	write(lun,97010) fromwhere, nstep, latndx(i), lonndx(i), k, ipair
!!$ 	write(lun,97020) 'drv old/oldbnd/new/del     ',   &
!!$ 		dryvol_t_old, dryvol_t_oldbnd, dryvol_t_new, dryvol_t_del
!!$ 	write(lun,97020) 'num old/oldbnd, dgnold/new ',   &
!!$ 		num_t_old, num_t_oldbnd, dgn_t_old, dgn_t_new
!!$ 	write(lun,97020) 'tailfr v_old/new, n_old/new',   &
!!$ 		tailfr_volold, tailfr_volnew, tailfr_numold, tailfr_numnew
!!$ 	dum = max(1.0e-10_r8,xferfrac_vol) / max(1.0e-10_r8,xferfrac_num)
!!$ 	dgn_xfer = dgn_t_new * dum**onethird
!!$ 	dum = max(1.0e-10_r8,(1.0_r8-xferfrac_vol)) /   &
!!$               max(1.0e-10_r8,(1.0_r8-xferfrac_num))
!!$ 	dgn_aftr = dgn_t_new * dum**onethird
!!$ 	write(lun,97020) 'xferfrac_v/n; dgn_xfer/aftr',   &
!!$ 		xferfrac_vol, xferfrac_num, dgn_xfer, dgn_aftr
!!$ !97010	format( / 'RENAME ', a, '  nx,lc,i,k,ip', i8, 4i4 )
!!$ 97010	format( / 'RENAME ', a, '  nx,lat,lon,k,ip', i8, 4i4 )
!!$ 97020	format( a, 6(1pe15.7) )
!!$ 	end if
!!$ 	end if
!   diagnostic output end   ------------------------------------------


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
8100	write(lunout,9050) ipair
	call endrun( 'modal_aero_rename_no_acc_crs_sub error' )
9050	format( / '*** subr. modal_aero_rename_no_acc_crs_sub ***' /   &
      	    4x, 'aerosol renaming not implemented for ipair =', i5 )

!EOC
	end subroutine modal_aero_rename_no_acc_crs_sub



!-------------------------------------------------------------------------
	subroutine modal_aero_rename_no_acc_crs_init
!
!   computes pointers for species transfer during aerosol renaming
!	(a2 --> a1 transfer)
!   transfers include number_a, number_c, mass_a, mass_c and
!	water_a
!

	implicit none

!   local variables
	integer :: ipair, iq, iqfrm, iqtoo
	integer :: lsfrma, lsfrmc, lstooa, lstooc, lunout
	integer :: mfrm, mtoo
	integer :: n1, n2, nspec
	integer :: nchfrma, nchfrmc, nchfrmskip, nchtooa, nchtooc, nchtooskip

	lunout = iulog
!
!   define "from mode" and "to mode" for each tail-xfer pairing
!	currently just a2-->a1
!
	n1 = modeptr_accum
	n2 = modeptr_aitken
	if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
	    npair_renamexf = 1
	    modefrm_renamexf(1) = n2
	    modetoo_renamexf(1) = n1
	else
	    npair_renamexf = 0
	    return
	end if

!
!   define species involved in each tail-xfer pairing
!	(include aerosol water)
!
aa_ipair: do ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)
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
		lsfrma = numptr_amode(mfrm)
		lstooa = numptr_amode(mtoo)
		lsfrmc = numptrcw_amode(mfrm)
		lstooc = numptrcw_amode(mtoo)
	    else if (iqfrm == 0) then
!   bypass transfer of aerosol water due to renaming
                cycle aa_iqfrm
!               lsfrma = lwaterptr_amode(mfrm)
!               lsfrmc = 0
!               lstooa = lwaterptr_amode(mtoo)
!               lstooc = 0
	    else
		lsfrma = lmassptr_amode(iqfrm,mfrm)
		lsfrmc = lmassptrcw_amode(iqfrm,mfrm)
		lstooa = 0
		lstooc = 0
	    end if


	    if ((lsfrma < 1) .or. (lsfrma > pcnst)) then
		write(lunout,9100) mfrm, iqfrm, lsfrma
		call endrun( 'modal_aero_rename_init error aa' )
	    end if
	    if ((lsfrmc < 1) .or. (lsfrmc > pcnst)) then
		write(lunout,9102) mfrm, iqfrm, lsfrmc
		call endrun( 'modal_aero_rename_init error bb' )
	    end if


	    if (iqfrm > 0) then
		nchfrma = len( trim( cnst_name(lsfrma) ) ) - nchfrmskip

! find "too" species having same lspectype_amode as the "frm" species
! AND same cnst_name (except for last 1/2/3 characters which are the mode index)
		do iqtoo = 1, nspec_amode(mtoo)
!		    if ( lspectype_amode(iqtoo,mtoo) .eq.   &
!			 lspectype_amode(iqfrm,mfrm) ) then
			lstooa = lmassptr_amode(iqtoo,mtoo)
			nchtooa = len( trim( cnst_name(lstooa) ) ) - nchtooskip
			if (cnst_name(lsfrma)(1:nchfrma) == cnst_name(lstooa)(1:nchtooa)) then
			! interstitial names match, so check cloudborne names too
			    nchfrmc = len( trim( cnst_name_cw(lsfrmc) ) ) - nchfrmskip
			    lstooc = lmassptrcw_amode(iqtoo,mtoo)
			    nchtooc = len( trim( cnst_name_cw(lstooc) ) ) - nchtooskip
			    if (cnst_name_cw(lsfrmc)(1:nchfrmc) /= &
			        cnst_name_cw(lstooc)(1:nchtooc)) lstooc = 0
			    exit
			else
			    lstooa = 0
			end if
!		    end if
		end do
	    end if ! (iqfrm > 0)

	    if ((lstooc < 1) .or. (lstooc > pcnst)) lstooc = 0
	    if ((lstooa < 1) .or. (lstooa > pcnst)) lstooa = 0
	    if (lstooa == 0) then
		write(lunout,9104) mfrm, iqfrm, lsfrma, iqtoo, lstooa
		call endrun( 'modal_aero_rename_init error cc' )
	    end if
	    if ((lstooc == 0) .and. (iqfrm /= 0)) then
		write(lunout,9104) mfrm, iqfrm, lsfrmc, iqtoo, lstooc
		call endrun( 'modal_aero_rename_init error dd' )
	    end if

	    nspec = nspec + 1
	    lspecfrma_renamexf(nspec,ipair) = lsfrma
	    lspectooa_renamexf(nspec,ipair) = lstooa
	    lspecfrmc_renamexf(nspec,ipair) = lsfrmc
	    lspectooc_renamexf(nspec,ipair) = lstooc
	end do aa_iqfrm

	nspecfrm_renamexf(ipair) = nspec
	end do aa_ipair

9100	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspecfrma out of range' /   &
      	'modefrm, ispecfrm, lspecfrma =', 3i6 / )
9102	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspecfrmc out of range' /   &
      	'modefrm, ispecfrm, lspecfrmc =', 3i6 / )
9104	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspectooa out of range' /   &
      	'modefrm, ispecfrm, lspecfrma, ispectoo, lspectooa =', 5i6 / )
9106	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspectooc out of range' /   &
      	'modefrm, ispecfrm, lspecfrmc, ispectoo, lspectooc =', 5i6 / )

!
!   output results
!
	if ( masterproc ) then

	write(lunout,9310)

	do 2900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)
	write(lunout,9320) ipair, mfrm, mtoo

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

2900	continue
	write(lunout,*)

	end if ! ( masterproc )

9310	format( / 'subr. modal_aero_rename_no_acc_crs_init' )
9320	format( 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350	format( 5x, 'no corresponding activated species' )

	return
	end subroutine modal_aero_rename_no_acc_crs_init

!----------------------------------------------------------------------
! code for troposphere and stratosphere
! -- allows accumulation to coarse mode exchange
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_rename_acc_crs_sub --- ...
!
! !INTERFACE:
	subroutine modal_aero_rename_acc_crs_sub(                       &
                        fromwhere,         lchnk,               &
                        ncol,              nstep,               &
                        loffset,           deltat,              &
                        pdel,              troplev,             &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx,           &
                        dqdt_rnpos                              )

! !USES:

   use physconst, only: gravit, mwdry
   use units, only: getunit
   use shr_spfn_mod, only: erfc => shr_spfn_erfc

   implicit none


! !PARAMETERS:
   character(len=*), intent(in) :: fromwhere    ! identifies which module
                                                ! is making the call
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)
   integer,  intent(in)    :: troplev(pcols)

   real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: q(ncol,pver,pcnstxx) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: qqcw(ncol,pver,pcnstxx) ! like q but for cloud-borne species

   real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the 
                              !     "fromwhere" continuous growth process 
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
   real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx)  
                              ! tendencies for "other" continuous growth process 
                              ! currently in cam3
                              !     dqdt is from gas (h2so4, nh3) condensation
                              !     dqdt_other is from aqchem and soa
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx)  
   logical,  intent(inout) :: dotendrn(pcnstxx) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(pcnstxx)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(ncol,pver) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
                              ! process-specific column tracer tendencies 
   real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
   real(r8), optional, intent(out) &
                           :: dqdt_rnpos(ncol,pver,pcnstxx)
                              ! the positive (production) part of the renaming tendency

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
   integer, save :: lun = -1  ! logical unit for diagnostics (6, or other
                              ! if a special diagnostics file is opened)

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

!   get logical unit (for output to dumpconv, deactivate the "lun = 6")
 	lun = iulog
	if (lun < 1) then
	   lun = getunit()
 	   open( unit=lun, file='dump.rename',   &
 			status='unknown', form='formatted' )
	end if


!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
!	if (npair_renamexf .eq. -123456789) then
!	    npair_renamexf = 0
!	    call modal_aero_rename_init
!	end if

!
!   check if any renaming pairs exist
!
	if (npair_renamexf .le. 0) return
! 	if (ncol .ne. -123456789) return
!	if (fromwhere .eq. 'aqchem') return


	deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
	onethird = 1.0_r8/3.0_r8
	xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

	if ( present( dqdt_rnpos ) ) then
	    l_dqdt_rnpos = .true.
	    dqdt_rnpos(:,:,:) = 0.0_r8
	else
	    l_dqdt_rnpos = .false.
	end if



!
!   loop over renaming pairs
!
mainloop1_ipair:  do ipair = 1, npair_renamexf

	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)

	flagaa_shrink = .false.
	if ((mfrm==modeptr_coarse) .and. (mtoo==modeptr_accum)) &
	    flagaa_shrink = .true.

!
!   compute aerosol dry-volume for the "from mode" of each renaming pair
!   also compute dry-volume change during the continuous growth process
!	using the incoming dqdt*deltat
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
!	in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
!	(old/new = before/after the continuous growth)
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


!!   diagnostic output start ----------------------------------------
!!	if (ldiag1 > 0) then
! 	icol_diag = -1
! 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) icol_diag = i
!!	if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
!! qak
! 	if (ldiag1 <= 0) then
! 	if ((i == 1) .and. (k == 1)) then
!! qak
! !	write(lun,97010) fromwhere, nstep, lchnk, i, k, ipair
!! 	write(lun,97010) fromwhere, nstep, latndx(i), lonndx(i), k, ipair
!! 	write(lun,97020) 'drv old/oldb/oldbnd/new/del    ',   &
!! 		dryvol_t_old, dryvol_t_oldb, dryvol_t_oldbnd, &
!! 		dryvol_t_new, dryvol_t_del
!! 	write(lun,97020) 'num old/oldbnd, dgnold/oldb/new',   &
!! 		num_t_old, num_t_oldbnd, dgn_t_old, dgn_t_oldb, dgn_t_new
!! 	write(lun,97020) 'tailfr v_old/new, n_old/new    ',   &
!! 		tailfr_volold, tailfr_volnew, tailfr_numold, tailfr_numnew
! 	tmpa = max(1.0e-10_r8,xferfrac_vol) / max(1.0e-10_r8,xferfrac_num)
! 	dgn_xfer = dgn_t_new * tmpa**onethird
! 	tmpa = max(1.0e-10_r8,(1.0_r8-xferfrac_vol)) /   &
!               max(1.0e-10_r8,(1.0_r8-xferfrac_num))
!! 	dgn_aftr = dgn_t_new * tmpa**onethird
!! 	write(lun,97020) 'xferfrac_v/n; dgn_xfer/aftr    ',   &
!! 		xferfrac_vol, xferfrac_num, dgn_xfer, dgn_aftr
! !97010	format( / 'RENAME ', a, '  nx,lc,i,k,ip', i8, 4i4 )
! 97010	format( / 'RENAME ', a, '  nx,lat,lon,k,ip', i8, 4i4 )
! 97020	format( a, 6(1pe15.7) )
! 	end if
! 	end if
!   diagnostic output end   ------------------------------------------


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
8100	write(lunout,9050) ipair
	call endrun( 'modal_aero_rename_acc_crs_sub error' )
9050	format( / '*** subr. modal_aero_rename_acc_crs_sub ***' /   &
      	    4x, 'aerosol renaming not implemented for ipair =', i5 )

!EOC
	end subroutine modal_aero_rename_acc_crs_sub



!-------------------------------------------------------------------------
! for modal aerosols in the troposphere and stratophere
! -- allows accumulation to coarse mode exchange
!-------------------------------------------------------------------------
	subroutine modal_aero_rename_acc_crs_init
!
!   computes pointers for species transfer during aerosol renaming
!	(a2 --> a1 transfer)
!   transfers include number_a, number_c, mass_a, mass_c and
!	water_a
!

	implicit none

!   local variables
	integer :: i, ipair, iq, iqfrm, iqtooa, iqtooc, itmpa
	integer :: l, lsfrma, lsfrmc, lstooa, lstooc, lunout
	integer :: mfrm, mtoo
	integer :: n1, n2, nspec
	integer :: nch_lfrm, nch_ltoo, nch_mfrmid, nch_mtooid

	real (r8) :: frelax

	lunout = iulog

!
!   define "from mode" and "to mode" for each tail-xfer pairing
!	using the values in ipair_select_renamexf(:)
!
	npair_renamexf = 0
	do ipair = 1, maxpair_renamexf
	    itmpa = ipair_select_renamexf(ipair)
	    if (itmpa == 0) then
		exit
	    else if (itmpa == 2001) then
		mfrm = modeptr_aitken
		mtoo = modeptr_accum
		igrow_shrink_renamexf(ipair) = 1
		ixferable_all_needed_renamexf(ipair) = 1
                strat_only_renamexf(ipair) = .false.
	    else if (itmpa == 1003) then
		mfrm = modeptr_accum
		mtoo = modeptr_coarse
		igrow_shrink_renamexf(ipair) = 1
		ixferable_all_needed_renamexf(ipair) = 0
                strat_only_renamexf(ipair) = .true.
	    else if (itmpa == 3001) then
		mfrm = modeptr_coarse
		mtoo = modeptr_accum
		igrow_shrink_renamexf(ipair) = -1
		ixferable_all_needed_renamexf(ipair) = 0
                strat_only_renamexf(ipair) = .true.
	    else
		write(lunout,'(/2a,3(1x,i12))') &
		    '*** subr. modal_aero_rename_acc_crs_init', &
		    'bad ipair_select_renamexf', ipair, itmpa
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if

	    do i = 1, ipair-1
		if (itmpa .eq. ipair_select_renamexf(i)) then
		    write(lunout,'(/2a/10(1x,i12))') &
			'*** subr. modal_aero_rename_acc_crs_init', &
			'duplicates in ipair_select_renamexf', &
			ipair_select_renamexf(1:ipair)
		    call endrun( 'modal_aero_rename_acc_crs_init error' )
		end if
	    end do

	    if ( (mfrm .ge. 1) .and. (mfrm .le. ntot_amode) .and. &
	         (mtoo .ge. 1) .and. (mtoo .le. ntot_amode) ) then
		npair_renamexf = ipair
		modefrm_renamexf(ipair) = mfrm
		modetoo_renamexf(ipair) = mtoo
	    else
		write(lunout,'(/2a,3(1x,i12))') &
		    '*** subr. modal_aero_rename_acc_crs_init', &
		    'bad mfrm or mtoo', ipair, mfrm, mtoo
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if
	end do ! ipair

	if (npair_renamexf .le. 0) then
	    write(lunout,'(/a/a,3(1x,i12))') &
		'*** subr. modal_aero_rename_acc_crs_init -- npair_renamexf = 0'
	    return
	end if


!
!   define species involved in each tail-xfer pairing
!	(include aerosol water)
!
	do 1900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)
	ixferable_all_renamexf(ipair) = 1

	if (mfrm < 10) then
	    nch_mfrmid = 1
	else if (mfrm < 100) then
	    nch_mfrmid = 2
	else
	    nch_mfrmid = 3
	end if
	if (mtoo < 10) then
	    nch_mtooid = 1
	else if (mtoo < 100) then
	    nch_mtooid = 2
	else
	    nch_mtooid = 3
	end if

	nspec = 0
	do 1490 iqfrm = -1, nspec_amode(mfrm)
	    if (iqfrm .eq. -1) then
		lsfrma = numptr_amode(mfrm)
		lstooa = numptr_amode(mtoo)
		lsfrmc = numptrcw_amode(mfrm)
		lstooc = numptrcw_amode(mtoo)
	    else if (iqfrm .eq. 0) then
!   bypass transfer of aerosol water due to renaming
                goto 1490
!               lsfrma = lwaterptr_amode(mfrm)
!               lsfrmc = 0
!               lstooa = lwaterptr_amode(mtoo)
!               lstooc = 0
	    else
		lsfrma = lmassptr_amode(iqfrm,mfrm)
		lsfrmc = lmassptrcw_amode(iqfrm,mfrm)
		lstooa = 0
		lstooc = 0
	    end if

	    if ((lsfrma .lt. 1) .or. (lsfrma .gt. pcnst)) then
		write(lunout,9100) ipair, mfrm, iqfrm, lsfrma
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if
	    if (iqfrm .le. 0) goto 1430

	    if ((lsfrmc .lt. 1) .or. (lsfrmc .gt. pcnst)) then
		write(lunout,9102) ipair, mfrm, iqfrm, lsfrmc
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if

! find "too" species having same name (except for mode number) as the "frm" species
	    nch_lfrm = len(trim(cnst_name(lsfrma))) - nch_mfrmid
	    iqtooa = -99
	    do iq = 1, nspec_amode(mtoo)
		l = lmassptr_amode(iq,mtoo)
		if ((l .lt. 1) .or. (l .gt. pcnst)) cycle
		nch_ltoo = len(trim(cnst_name(l))) - nch_mtooid
		if ( cnst_name(lsfrma)(1:nch_lfrm) == &
		     cnst_name(l     )(1:nch_ltoo) ) then
		    lstooa = l
		    iqtooa = iq
		    exit
		end if
	    end do

	    nch_lfrm = len(trim(cnst_name_cw(lsfrmc))) - nch_mfrmid
	    iqtooc = -99
	    do iq = 1, nspec_amode(mtoo)
		l = lmassptrcw_amode(iq,mtoo)
		if ((l .lt. 1) .or. (l .gt. pcnst)) cycle
		nch_ltoo = len(trim(cnst_name_cw(l))) - nch_mtooid
		if ( cnst_name_cw(lsfrmc)(1:nch_lfrm) == &
		     cnst_name_cw(l     )(1:nch_ltoo) ) then
		    lstooc = l
		    iqtooc = iq
		    exit
		end if
	    end do

1430	    if ((lstooc .lt. 1) .or. (lstooc .gt. pcnst)) lstooc = 0
	    if ((lstooa .lt. 1) .or. (lstooa .gt. pcnst)) lstooa = 0

	    if ((lstooa .eq. 0) .or. (lstooc .eq. 0)) then
		if ( ( masterproc                                  ) .or. &
		     ( (lstooa .ne. 0) .or. (lstooc .ne. 0)        ) .or. &
		     ( ixferable_all_needed_renamexf(ipair) .gt. 0 ) ) then
		    if (lstooa .eq. 0) &
			write(lunout,9104) trim(cnst_name(lsfrma)), &
			    ipair, mfrm, iqfrm, lsfrma, iqtooa, lstooa
		    if (lstooc .eq. 0) &
			write(lunout,9106) trim(cnst_name_cw(lsfrmc)), &
			    ipair, mfrm, iqfrm, lsfrmc, iqtooc, lstooc
		end if
		if ((lstooa .ne. 0) .or. (lstooc .ne. 0)) then
		    write(lunout,9108)
		    call endrun( 'modal_aero_rename_acc_crs_init error' )
		end if
		if (ixferable_all_needed_renamexf(ipair) .gt. 0) then
		    write(lunout,9109)
		    call endrun( 'modal_aero_rename_acc_crs_init error' )
		end if
		ixferable_all_renamexf(ipair) = 0
		if (iqfrm .gt. 0) then
		    ixferable_a_renamexf(iqfrm,ipair) = 0
		    ixferable_c_renamexf(iqfrm,ipair) = 0
		end if
	    else
		nspec = nspec + 1
		lspecfrma_renamexf(nspec,ipair) = lsfrma
		lspectooa_renamexf(nspec,ipair) = lstooa
		lspecfrmc_renamexf(nspec,ipair) = lsfrmc
		lspectooc_renamexf(nspec,ipair) = lstooc
		if (iqfrm .gt. 0) then
		    ixferable_a_renamexf(iqfrm,ipair) = 1
		    ixferable_c_renamexf(iqfrm,ipair) = 1
		end if
	    end if
1490	continue

	nspecfrm_renamexf(ipair) = nspec
1900	continue

9100	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspecfrma out of range' /   &
      	'ipair, modefrm, ispecfrm, lspecfrma =', 4i6 )
9102	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspecfrmc out of range' /   &
      	'ipair, modefrm, ispecfrm, lspecfrmc =', 4i6 )
9104	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspectooa out of range for', 2x, a /   &
      	'ipair, modefrm, ispecfrm, lspecfrma, ispectoo, lspectooa =', 6i6 )
9106	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspectooc out of range for', 2x, a /   &
      	'ipair, modefrm, ispecfrm, lspecfrmc, ispectoo, lspectooc =', 6i6 )
9108	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'only one of lspectooa and lspectooc is out of range' )
9109	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'all species must be xferable for this pair' )


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

	    dp_cut(ipair) = sqrt(   &
		dgnum_amode(mfrm)*exp(1.5_r8*(alnsg_amode(mfrm)**2)) *   &
		dgnum_amode(mtoo)*exp(1.5_r8*(alnsg_amode(mtoo)**2)) )
	    dp_xferall_thresh(ipair) = dgnum_amode(mtoo)
	    dp_xfernone_threshaa(ipair) = dgnum_amode(mfrm)

	    if ((mfrm == modeptr_accum) .and. (mtoo == modeptr_coarse)) then
		dp_cut(ipair)               = 4.4e-7_r8 
		dp_xfernone_threshaa(ipair) = 1.6e-7_r8 
		dp_xferall_thresh(ipair)    = 4.7e-7_r8 
	    else if ((mfrm == modeptr_coarse) .and. (mtoo == modeptr_accum)) then
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

2900	continue
	write(lunout,*)

	end if ! ( masterproc )

9310	format( / 'subr. modal_aero_rename_acc_crs_init' )
9320	format( / 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3, &
	        5x, 'igrow_shrink', i3, 5x, 'ixferable_all', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350	format( 5x, 'no corresponding activated species' )


	return
	end subroutine modal_aero_rename_acc_crs_init

!----------------------------------------------------------------------

   end module modal_aero_rename
