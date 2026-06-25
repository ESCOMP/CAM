! CAM wrapper for modal_aero_rename.
! Owns the resolved renaming-pair tables (consumed by other CAM aerosol
! wrappers), performs the cnst_name-based pair resolution, and hands the
! tables + mode metadata to the portable modal_aero_rename_init.
!----------------------------------------------------------------------
  module modal_aero_rename_cam

! !USES:
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use cam_abortutils,  only: endrun
  use cam_logfile,     only: iulog
  use mo_constants,    only: pi
  use constituents,    only: pcnst, cnst_name
  use spmd_utils,      only: masterproc
  use modal_aero_data, only: maxspec_renamexf=>nspec_max, ntot_amode
  use modal_aero_data, only: alnsg_amode, voltonumblo_amode, voltonumbhi_amode, dgnum_amode, nspec_amode
  use modal_aero_data, only: lmassptr_amode, lmassptrcw_amode
  use modal_aero_data, only: numptr_amode, numptrcw_amode, modeptr_coarse, modeptr_accum
  use modal_aero_data, only: modeptr_stracoar
  use modal_aero_data, only: dgnumhi_amode, dgnumlo_amode, cnst_name_cw, modeptr_aitken
  use radiative_aerosol,only: rad_aer_get_mode_idx
  use modal_aero_rename, only: maxpair_renamexf

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public :: modal_aero_rename_cam_init

! !PUBLIC DATA MEMBERS:
! Resolved renaming-pair tables (host constituent-index space).  These are
! consumed by modal_aero_calcsize_cam and modal_aero_gasaerexch_cam.
  integer, public :: npair_renamexf = -123456789
  integer, protected, public :: modefrm_renamexf(maxpair_renamexf)
  integer, protected, public :: modetoo_renamexf(maxpair_renamexf)
  integer, protected, public :: nspecfrm_renamexf(maxpair_renamexf)

  integer, allocatable, protected, public :: lspecfrma_renamexf(:,:)
  integer, allocatable, protected, public :: lspecfrmc_renamexf(:,:)
  integer, allocatable, protected, public :: lspectooa_renamexf(:,:)
  integer, allocatable, protected, public :: lspectooc_renamexf(:,:)

! ipair_select_renamexf defines the mode_from and mode_too for each renaming pair
! 2001 = aitken --> accum
! 1003 = accum  --> coarse
! 3001 = coarse --> accum
! 1005 = accum  --> stracoar
! 5001 = stracoar --> accum
  integer :: ipair_select_renamexf(maxpair_renamexf)

! Renaming-pair flags resolved here and consumed by the portable science: handed
! to modal_aero_rename_init, and passed by aero_model to modal_aero_rename_run.
  integer, protected, public :: igrow_shrink_renamexf(maxpair_renamexf)
  integer, protected, public :: ixferable_all_renamexf(maxpair_renamexf)
  integer :: ixferable_all_needed_renamexf(maxpair_renamexf)
  integer, allocatable, protected, public :: ixferable_a_renamexf(:,:)
  integer, allocatable, protected, public :: ixferable_c_renamexf(:,:)

  logical, protected, public :: strat_only_renamexf(maxpair_renamexf)
! strat_only_renamexf - when true for a particular renaming pair, renaming is only
!                       done in stratosphere (when k < troplev(icol) )

  logical :: modal_accum_coarse_exch = .false.

!----------------------------------------------------------------------
contains

  !------------------------------------------------------------------
  ! Resolve the renaming-pair tables from CAM constituent metadata, then hand
  ! them (with the mode metadata) to the portable modal_aero_rename_init.
  !------------------------------------------------------------------
  subroutine modal_aero_rename_cam_init(modal_accum_coarse_exch_in)
    use modal_aero_rename, only: modal_aero_rename_init

    logical, optional, intent(in) :: modal_accum_coarse_exch_in

    character(len=512) :: errmsg
    integer            :: errflg

    ! ipair_select_renamexf defines the mode_from and mode_too for each renaming pair
    ! 2001 = aitken --> accum
    ! 1003 = accum  --> coarse
    ! 3001 = coarse --> accum
    ! 1005 = accum  --> stracoar
    ! 5001 = stracoar --> accum
    if( rad_aer_get_mode_idx(0,'coarse_strat') > 0 ) then
       ipair_select_renamexf(1:maxpair_renamexf) = (/ 2001, 1005, 5001 /)
    else
       ipair_select_renamexf(1:maxpair_renamexf) = (/ 2001, 1003, 3001 /)
    endif

    allocate( lspecfrma_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( lspecfrmc_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( lspectooa_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( lspectooc_renamexf(maxspec_renamexf,maxpair_renamexf) )

    allocate( ixferable_a_renamexf(maxspec_renamexf,maxpair_renamexf) )
    allocate( ixferable_c_renamexf(maxspec_renamexf,maxpair_renamexf) )

    ! Default-initialize the accum-coarse-exchange flags: the no_acc_crs
    ! resolution path leaves them unset, and they are handed (unused) to the
    ! portable code on that path.
    igrow_shrink_renamexf(:)  = 0
    ixferable_all_renamexf(:) = 0
    ixferable_a_renamexf(:,:) = 0
    ixferable_c_renamexf(:,:) = 0
    strat_only_renamexf(:)    = .false.

    if (present(modal_accum_coarse_exch_in)) then
       modal_accum_coarse_exch = modal_accum_coarse_exch_in
    endif

    if (modal_accum_coarse_exch) then
       call modal_aero_rename_acc_crs_init()
    else
       call modal_aero_rename_no_acc_crs_init()
    endif

    ! Precompute rename's own accum-coarse-exchange physics coefficients. Only
    ! the metadata the precompute + one-time log needs is handed over; the shared
    ! tables/metadata are passed per-call to modal_aero_rename_run instead.
    call modal_aero_rename_init(                                     &
       modal_accum_coarse_exch = modal_accum_coarse_exch,            &
       ntot_amode            = ntot_amode,                           &
       alnsg_amode           = alnsg_amode,                          &
       dgnum_amode           = dgnum_amode,                          &
       dgnumhi_amode         = dgnumhi_amode,                        &
       dgnumlo_amode         = dgnumlo_amode,                        &
       voltonumblo_amode     = voltonumblo_amode,                    &
       voltonumbhi_amode     = voltonumbhi_amode,                    &
       modeptr_accum         = modeptr_accum,                        &
       modeptr_coarse        = modeptr_coarse,                       &
       modeptr_stracoar      = modeptr_stracoar,                     &
       npair_renamexf        = npair_renamexf,                       &
       modefrm_renamexf      = modefrm_renamexf,                     &
       modetoo_renamexf      = modetoo_renamexf,                     &
       nspecfrm_renamexf     = nspecfrm_renamexf,                    &
       lspecfrma_renamexf    = lspecfrma_renamexf,                   &
       lspecfrmc_renamexf    = lspecfrmc_renamexf,                   &
       lspectooa_renamexf    = lspectooa_renamexf,                   &
       lspectooc_renamexf    = lspectooc_renamexf,                   &
       igrow_shrink_renamexf  = igrow_shrink_renamexf,               &
       ixferable_all_renamexf = ixferable_all_renamexf,              &
       cnst_name_in          = cnst_name,                            &
       cnst_name_cw_in       = cnst_name_cw,                         &
       pi                    = pi,                                   &
       amRoot                = masterproc,                           &
       iulog_in              = iulog,                                &
       errmsg                = errmsg,                               &
       errflg                = errflg                                )

    if (errflg /= 0) then
       call endrun('modal_aero_rename_cam_init: '//trim(errmsg))
    end if

  end subroutine modal_aero_rename_cam_init

!----------------------------------------------------------------------
! private methods -- renaming-pair resolution (verbatim from the original
! modal_aero_rename init routines; cnst_name matching is host-specific)
!----------------------------------------------------------------------

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
           else if (itmpa == 2001) then  !both mam4 and mam5
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
           else if (itmpa == 1005) then
              mfrm = modeptr_accum
              mtoo = modeptr_stracoar
              igrow_shrink_renamexf(ipair) = 1
              ixferable_all_needed_renamexf(ipair) = 0
              strat_only_renamexf(ipair) = .true.
           else if (itmpa == 3001) then
              mfrm = modeptr_coarse
              mtoo = modeptr_accum
              igrow_shrink_renamexf(ipair) = -1
              ixferable_all_needed_renamexf(ipair) = 0
              strat_only_renamexf(ipair) = .true.
           else if (itmpa == 5001) then
              mfrm = modeptr_stracoar
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
	end subroutine modal_aero_rename_acc_crs_init

  end module modal_aero_rename_cam
