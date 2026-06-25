! CAM wrapper for modal_aero_coag.
! Owns the pair_option_acoag build-time selection, resolves the
! coagulation-pair tables from CAM constituent metadata (with the pcage
! aging tables from modal_aero_gasaerexch), hands them (with the mode
! metadata and host physical constants) to the portable
! modal_aero_coag_init, and registers history fields.
!----------------------------------------------------------------------
  module modal_aero_coag_cam

! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public :: modal_aero_coag_cam_init

#if ( defined MODAL_AERO_7MODE || defined MODAL_AERO_4MODE || defined MODAL_AERO_5MODE)
  integer, parameter :: pair_option_acoag = 3
#elif ( defined MODAL_AERO_3MODE )
  integer, parameter :: pair_option_acoag = 1
#endif
! specifies pairs of modes for which coagulation is calculated
!   1 -- [aitken-->accum]
!   2 -- [aitken-->accum], and [pcarbon-->accum]
!   3 -- [aitken-->accum], [pcarbon-->accum],
!        and [aitken-->pcarbon--(aging)-->accum]
! other -- do no coag

!----------------------------------------------------------------------
contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
	subroutine modal_aero_coag_cam_init
!
!   computes pointers for species transfer during coagulation
!   and hands them to the portable modal_aero_coag_init
!
	use modal_aero_coag, only:  modal_aero_coag_init, maxpair_acoag
	use modal_aero_data
	use modal_aero_gasaerexch, only:  &
		modefrm_pcage, nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage, &
                soa_equivso4_factor

	use cam_abortutils,  only: endrun
	use cam_history,     only: addfld, add_default, fieldname_len, horiz_only
	use constituents,    only: pcnst, cnst_name
	use physconst,       only: r_universal, pstd, tmelt, boltz
	use spmd_utils,      only: masterproc
        use phys_control,    only: phys_getopts

	implicit none

!   local variables
	integer :: ipair, iq, iqfrm, iqfrm_aa, iqtoo, iqtoo_aa
        integer :: jsoa
	integer :: l, l1, l2, lsfrm, lstoo, lunout
	integer :: m, mait, mpca, mfrm, mtoo, mtef
	integer :: nchfrm, nchfrmskip, nchtoo, nchtooskip, nspec

!   resolved coagulation-pair tables (host constituent-index space),
!   handed to the portable modal_aero_coag_init at the end
	integer :: maxspec_acoag
	integer :: npair_acoag
	integer :: modefrm_acoag(maxpair_acoag)
	integer :: modetoo_acoag(maxpair_acoag)
	integer :: modetooeff_acoag(maxpair_acoag)
	integer :: nspecfrm_acoag(maxpair_acoag)
	integer, allocatable :: lspecfrm_acoag(:,:)
	integer, allocatable :: lspectoo_acoag(:,:)
	integer :: ip_aitacc, ip_aitpca, ip_pcaacc
	real(r8), allocatable :: fac_m2v_aitage(:), fac_m2v_pcarbon(:)

	character(len=fieldname_len)   :: tmpname
	character(len=fieldname_len+3) :: fieldname
	character(128)                 :: long_name
	character(8)                   :: unit

	logical :: dotend(pcnst)
        logical :: history_aerosol      ! Output the MAM aerosol tendencies

	character(len=200) :: msg

	character(len=512) :: errmsg
	integer            :: errflg

        !-----------------------------------------------------------------------
        call phys_getopts( history_aerosol_out        = history_aerosol   )

        lunout = 6

        maxspec_acoag = nspec_max
        allocate( lspecfrm_acoag(maxspec_acoag,maxpair_acoag) )
        allocate( lspectoo_acoag(maxspec_acoag,maxpair_acoag) )
        allocate( fac_m2v_aitage(nspec_max), fac_m2v_pcarbon(nspec_max) )

!   default-initialize the tables so the unused slots are well-defined when
!   the whole arrays are handed to the portable modal_aero_coag_init
        modefrm_acoag(:)    = 0
        modetoo_acoag(:)    = 0
        modetooeff_acoag(:) = 0
        nspecfrm_acoag(:)   = 0
        lspecfrm_acoag(:,:) = 0
        lspectoo_acoag(:,:) = 0

!
!   define "from mode" and "to mode" for each coagulation pairing
!	currently just a2-->a1 coagulation
!
	if (pair_option_acoag == 1) then
	    npair_acoag = 1
	    modefrm_acoag(1) = modeptr_aitken
	    modetoo_acoag(1) = modeptr_accum
	    modetooeff_acoag(1) = modeptr_accum
	else if (pair_option_acoag == 2) then
	    npair_acoag = 2
	    modefrm_acoag(1) = modeptr_aitken
	    modetoo_acoag(1) = modeptr_accum
	    modetooeff_acoag(1) = modeptr_accum
	    modefrm_acoag(2) = modeptr_pcarbon
	    modetoo_acoag(2) = modeptr_accum
	    modetooeff_acoag(2) = modeptr_accum
	else if (pair_option_acoag == 3) then
	    npair_acoag = 3
	    modefrm_acoag(1) = modeptr_aitken
	    modetoo_acoag(1) = modeptr_accum
	    modetooeff_acoag(1) = modeptr_accum
	    modefrm_acoag(2) = modeptr_pcarbon
	    modetoo_acoag(2) = modeptr_accum
	    modetooeff_acoag(2) = modeptr_accum
	    modefrm_acoag(3) = modeptr_aitken
	    modetoo_acoag(3) = modeptr_pcarbon
	    modetooeff_acoag(3) = modeptr_accum
	    if (modefrm_pcage <= 0) then
		write(*,*) '*** modal_aero_coag_init error'
		write(*,*) '    pair_option_acoag, modefrm_pcage mismatch'
		write(*,*) '    pair_option_acoag, modefrm_pcage =', &
		    pair_option_acoag, modefrm_pcage
		call endrun( 'modal_aero_coag_init error' )
	    end if
	else
	    npair_acoag = 0
	    return
	end if

!
!   define species involved in each coagulation pairing
!	(include aerosol water)
!
aa_ipair: do ipair = 1, npair_acoag

	mfrm = modefrm_acoag(ipair)
	mtoo = modetoo_acoag(ipair)
	mtef = modetooeff_acoag(ipair)
	if ( (mfrm < 1) .or. (mfrm > ntot_amode) .or.   &
	     (mtoo < 1) .or. (mtoo > ntot_amode) .or.   &
	     (mtef < 1) .or. (mtef > ntot_amode) ) then
	    write(*,*) '*** modal_aero_coag_init error'
	    write(*,*) '    ipair, ntot_amode =', ipair, ntot_amode
	    write(*,*) '    mfrm, mtoo, mtef  =', mfrm, mtoo, mtef
	    call endrun( 'modal_aero_coag_init error' )
	end if


	mtoo = mtef   ! effective modetoo
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
aa_iqfrm: do iqfrm = 1, nspec_amode(mfrm)
	    lsfrm = lmassptr_amode(iqfrm,mfrm)
	    if ((lsfrm .lt. 1) .or. (lsfrm .gt. pcnst)) cycle aa_iqfrm
	    nchfrm = len( trim( cnst_name(lsfrm) ) ) - nchfrmskip
! find "too" species having same lspectype_amode as the "frm" species
! AND same cnst_name (except for last 1/2/3 characters which are the mode index)
	    do iqtoo = 1, nspec_amode(mtoo)
              lstoo = lmassptr_amode(iqtoo,mtoo)
              nchtoo = len( trim( cnst_name(lstoo) ) ) - nchtooskip
              if (cnst_name(lsfrm)(1:nchfrm) == cnst_name(lstoo)(1:nchtoo)) then
                 exit
              else
                 lstoo = 0
              end if
	    end do

	    if ((lstoo < 1) .or. (lstoo > pcnst)) lstoo = 0
	    nspec = nspec + 1
	    lspecfrm_acoag(nspec,ipair) = lsfrm
	    lspectoo_acoag(nspec,ipair) = lstoo
	end do aa_iqfrm

!       lsfrm = lwaterptr_amode(mfrm)
!       if ((lsfrm .ge. 1) .and. (lsfrm .le. pcnst)) then
!           lstoo = lwaterptr_amode(mtoo)
!           if ((lstoo .lt. 1) .or. (lstoo .gt. pcnst)) lstoo = 0
!           nspec = nspec + 1
!           lspecfrm_acoag(nspec,ipair) = lsfrm
!           lspectoo_acoag(nspec,ipair) = lstoo
!       end if

	nspecfrm_acoag(ipair) = nspec
	end do aa_ipair

!
!   output results
!
	if ( masterproc ) then

	write(lunout,9310)

	do ipair = 1, npair_acoag
	  mfrm = modefrm_acoag(ipair)
	  mtoo = modetoo_acoag(ipair)
	  mtef = modetooeff_acoag(ipair)
	  write(lunout,9320) ipair, mfrm, mtoo, mtef

	  do iq = 1, nspecfrm_acoag(ipair)
	    lsfrm = lspecfrm_acoag(iq,ipair)
	    lstoo = lspectoo_acoag(iq,ipair)
	    if (lstoo .gt. 0) then
		write(lunout,9330) lsfrm, cnst_name(lsfrm),   &
      			lstoo, cnst_name(lstoo)
	    else
		write(lunout,9340) lsfrm, cnst_name(lsfrm)
	    end if
	  end do

	end do ! ipair = ...
	write(lunout,*)

	end if ! ( masterproc ) 

9310	format( / 'subr. modal_aero_coag_init' )
9320	format( 'pair', i3, 5x, 'mode', i3, &
		' ---> mode', i3, '   eff', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )

!   set following variables that are used in modal_aero_coag_subr
!
	fac_m2v_aitage(:) = 0.0_r8
	fac_m2v_pcarbon(:) = 0.0_r8
	if (pair_option_acoag == 3) then
!   following ipair definitions MUST BE CONSISTENT with
!   the coding in modal_aero_coag_init for pair_option_acoag == 3
	    ip_aitacc = 1
	    ip_pcaacc = 2
	    ip_aitpca = 3

	    mait = modeptr_aitken
	    mpca = modeptr_pcarbon

	    ipair = ip_aitpca
	    do iq = 1, nspecfrm_acoag(ipair)
		lsfrm = lspecfrm_acoag(iq,ipair)
		l2 = -1
		do l1 = 1, nspec_amode(mait)
		   if (lmassptr_amode(l1,mait) == lsfrm) then
                      l2 = l1
			exit
		   end if
		end do
		if (l2 <= 0) then
		    write( msg, '(a,5(1x,i12))' ) &
			'modal_aero_coag_init error a001 for ipair, iq, lsfrm', &
			ipair, iq, lsfrm
		    call endrun( msg )
		end if
		if (lsfrm == lptr_so4_a_amode(mait)) then
!		    fac_m2v_aitage(iq) = specmw_amode(l2) / specdens_amode(l2) 
                   fac_m2v_aitage(iq) = specmw_amode(l1,mait) / specdens_amode(l1,mait)
		else if (lsfrm == lptr_nh4_a_amode(mait)) then
!                   fac_m2v_aitage(iq) = specmw_amode(l2) / specdens_amode(l2)
                    fac_m2v_aitage(iq) = specmw_amode(l1,mait) / specdens_amode(l1,mait)
		else 
		    do jsoa = 1, nsoa
			if (lsfrm == lptr2_soa_a_amode(mait,jsoa)) then
			    fac_m2v_aitage(iq) = soa_equivso4_factor(jsoa)*   &
                                 !(specmw_amode(l2) / specdens_amode(l2))
                                 (specmw_amode(l1,mait) / specdens_amode(l1,mait))
			end if
!   for soa, the soa_equivso4_factor converts the soa volume into an
!	so4(+nh4) volume that has same hygroscopicity contribution as soa
!   this allows aging calculations to be done in terms of the amount
!	of (equivalent) so4(+nh4) in the shell
!   (see modal_aero_gasaerexch)
		    end do
		end if
	    end do
	    
	    do l = 1, nspec_amode(mpca)
!B		l2 = lspectype_amode(l,mpca)
!   fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!		[m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
!		fac_m2v_pcarbon(l) = specmw_amode(l2) / specdens_amode(l2)
                fac_m2v_pcarbon(l) = specmw_amode(l,mpca) / specdens_amode(l,mpca)
	    end do

	else
	    ip_aitacc = -999888777
	    ip_pcaacc = -999888777
	    ip_aitpca = -999888777
	end if

!   hand the resolved tables, mode metadata, and host physical constants
!   to the portable scheme
	call modal_aero_coag_init(                            &
	   pair_option_acoag_in = pair_option_acoag,          &
	   npair_acoag_in       = npair_acoag,                &
	   modefrm_acoag_in     = modefrm_acoag,              &
	   modetoo_acoag_in     = modetoo_acoag,              &
	   modetooeff_acoag_in  = modetooeff_acoag,           &
	   nspecfrm_acoag_in    = nspecfrm_acoag,             &
	   lspecfrm_acoag_in    = lspecfrm_acoag,             &
	   lspectoo_acoag_in    = lspectoo_acoag,             &
	   ip_aitacc_in         = ip_aitacc,                  &
	   ip_aitpca_in         = ip_aitpca,                  &
	   ip_pcaacc_in         = ip_pcaacc,                  &
	   fac_m2v_aitage_in    = fac_m2v_aitage,             &
	   fac_m2v_pcarbon_in   = fac_m2v_pcarbon,            &
	   nspec_max_in         = nspec_max,                  &
	   ntot_amode_in        = ntot_amode,                 &
	   modeptr_accum_in     = modeptr_accum,              &
	   modeptr_aitken_in    = modeptr_aitken,             &
	   modeptr_pcarbon_in   = modeptr_pcarbon,            &
	   numptr_amode_in      = numptr_amode,               &
	   mprognum_amode_in    = mprognum_amode,             &
	   nspec_amode_in       = nspec_amode,                &
	   lmassptr_amode_in    = lmassptr_amode,             &
	   alnsg_amode_in       = alnsg_amode,                &
	   sigmag_amode_in      = sigmag_amode,               &
	   r_universal_in       = r_universal,                &
	   pstd_in              = pstd,                       &
	   tmelt_in             = tmelt,                      &
	   boltz_in             = boltz,                      &
	   errmsg               = errmsg,                     &
	   errflg               = errflg                      )

	if (errflg /= 0) then
	   call endrun('modal_aero_coag_cam_init: '//trim(errmsg))
	end if

!
!   create history file column-tendency fields
!
	dotend(:) = .false.
	do ipair = 1, npair_acoag
	  do iq = 1, nspecfrm_acoag(ipair)
	    l = lspecfrm_acoag(iq,ipair)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	    l = lspectoo_acoag(iq,ipair)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	  end do

	  m = modefrm_acoag(ipair)
	  if ((m > 0) .and. (m <= ntot_amode)) then
	    l = numptr_amode(m)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	  end if
	  m = modetoo_acoag(ipair)
	  if ((m > 0) .and. (m <= ntot_amode)) then
	    l = numptr_amode(m)
	    if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
	  end if
	end do ! ipair = ...

	if (pair_option_acoag == 3) then
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

	do l = 1, pcnst
	    if ( .not. dotend(l) ) cycle
	    tmpname = cnst_name(l)
	    unit = 'kg/m2/s'
	    do m = 1, ntot_amode
	        if (l == numptr_amode(m)) unit = '#/m2/s'
	    end do
	    fieldname = trim(tmpname) // '_sfcoag1'
	    long_name = trim(tmpname) // ' modal_aero coagulation column tendency'
	    call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
	    endif
	    if ( masterproc ) write(*,'(3(a,2x))') &
		'modal_aero_coag_init addfld', fieldname, unit
	end do ! l = ...


	return
	end subroutine modal_aero_coag_cam_init

  end module modal_aero_coag_cam
