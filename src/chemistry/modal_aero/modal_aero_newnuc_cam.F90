! CAM wrapper for modal_aero_newnuc.
! Resolves CAM-specific species indices and aitken-mode so4/nh4 properties,
! hands them (with the host physical constants) to the portable
! modal_aero_newnuc_init, and registers history fields.
!----------------------------------------------------------------------
module modal_aero_newnuc_cam

! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:
  public :: modal_aero_newnuc_cam_init

!----------------------------------------------------------------------
contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine modal_aero_newnuc_cam_init

!-----------------------------------------------------------------------
!
! Purpose:
!    resolve the species indices and aitken-mode so4/nh4 properties
!       needed by modal_aero_newnuc and hand them to the portable
!       modal_aero_newnuc_init
!    create history fields for column tendencies associated with
!       modal_aero_newnuc
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

use modal_aero_data
use modal_aero_newnuc, only:  modal_aero_newnuc_init

use cam_abortutils,   only:  endrun
use cam_history,      only:  addfld, add_default, fieldname_len, horiz_only
use constituents,     only:  pcnst, cnst_get_ind, cnst_name
use mo_constants,     only:  pi, rgas, avogadro
use physconst,        only:  mwso4, mwnh4, r_universal
use spmd_utils,       only:  masterproc
use phys_control,     only: phys_getopts


implicit none

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
   integer  :: l_h2so4, l_nh3
   integer  :: lnumait, lnh4ait, lso4ait
   integer  :: l, l1, l2
   integer  :: m, mait

   real(r8) :: mw_so4a_host, mw_nh4a_host
   real(r8) :: dens_so4a_host

   character(len=fieldname_len)   :: tmpname
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit

   logical                        :: dotend(pcnst)
   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies

   character(len=512)             :: errmsg
   integer                        :: errflg

   !-----------------------------------------------------------------------

        call phys_getopts( history_aerosol_out        = history_aerosol   )


!   set these indices
!   skip if no h2so4 species (the portable module keeps its bypass
!      defaults when modal_aero_newnuc_init is not called)
!   skip if no aitken mode so4 or num species
	call cnst_get_ind( 'H2SO4', l_h2so4, .false. )
	call cnst_get_ind( 'NH3', l_nh3, .false. )

	mait = modeptr_aitken
	if (mait > 0) then
	    lnumait = numptr_amode(mait)
	    lso4ait = lptr_so4_a_amode(mait)
	    lnh4ait = lptr_nh4_a_amode(mait)
	end if
	if ((l_h2so4  <= 0) .or. (l_h2so4 > pcnst)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- l_h2so4 <= 0'
	    return
	else if ((lso4ait <= 0) .or. (lso4ait > pcnst)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- lso4ait <= 0'
	    return
	else if ((lnumait <= 0) .or. (lnumait > pcnst)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- lnumait <= 0'
	    return
	else if ((mait <= 0) .or. (mait > ntot_amode)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- modeptr_aitken <= 0'
	    return
	end if

!   set these constants
!      mw_so4a_host is molec-wght of sulfate aerosol in host code
!         96 when nh3/nh4 are simulated
!         something else when nh3/nh4 are not simulated
	l = lptr_so4_a_amode(mait) ; l2 = -1
        if (l <= 0) call endrun( 'modal_aero_newnuch_init error a001 finding aitken so4' )
        do l1 = 1, nspec_amode(mait)
           if (lmassptr_amode(l1,mait) == l) then
              l2 = l1
              mw_so4a_host = specmw_amode(l1,mait)
              dens_so4a_host = specdens_amode(l1,mait)
           end if
        end do
        if (l2 <= 0) call endrun( 'modal_aero_newnuch_init error a002 finding aitken so4' )

        l = lptr_nh4_a_amode(mait) ; l2 = -1
        if (l > 0) then
           do l1 = 1, nspec_amode(mait)
              if (lmassptr_amode(l1,mait) == l) then
                 l2 = l1
                 mw_nh4a_host = specmw_amode(l1,mait)
              end if
           end do
           if (l2 <= 0) call endrun( 'modal_aero_newnuch_init error a002 finding aitken nh4' )
        else
           mw_nh4a_host = mw_so4a_host
        end if

!   hand the resolved indices, aitken-mode properties, and host physical
!   constants to the portable scheme
	call modal_aero_newnuc_init(                          &
	   l_h2so4_in         = l_h2so4,                      &
	   l_nh3_in           = l_nh3,                        &
	   lnumait_in         = lnumait,                      &
	   lnh4ait_in         = lnh4ait,                      &
	   lso4ait_in         = lso4ait,                      &
	   mw_so4a_host_in    = mw_so4a_host,                 &
	   mw_nh4a_host_in    = mw_nh4a_host,                 &
	   dens_so4a_host_in  = dens_so4a_host,               &
	   dgnum_aitken_in    = dgnum_amode(mait),            &
	   dgnumhi_aitken_in  = dgnumhi_amode(mait),          &
	   dgnumlo_aitken_in  = dgnumlo_amode(mait),          &
	   pi_in              = pi,                           &
	   rgas_in            = rgas,                         &
	   avogad_in          = avogadro,                     &
	   mw_so4a_in         = mwso4,                        &
	   mw_nh4a_in         = mwnh4,                        &
	   r_universal_in     = r_universal,                  &
	   errmsg             = errmsg,                       &
	   errflg             = errflg                        )

	if (errflg /= 0) then
	   call endrun('modal_aero_newnuc_cam_init: '//trim(errmsg))
	end if

!
!   create history file column-tendency fields
!
	dotend(:) = .false.
	dotend(lnumait) = .true.
	dotend(lso4ait) = .true.
	dotend(l_h2so4) = .true.
	if ((l_nh3   > 0) .and. (l_nh3   <= pcnst) .and. &
	    (lnh4ait > 0) .and. (lnh4ait <= pcnst)) then
	    dotend(lnh4ait) = .true.
	    dotend(l_nh3) = .true.
	end if

	do l = 1, pcnst
	    if ( .not. dotend(l) ) cycle
	    tmpname = cnst_name(l)
	    unit = 'kg/m2/s'
	    do m = 1, ntot_amode
	        if (l == numptr_amode(m)) unit = '#/m2/s'
	    end do
	    fieldname = trim(tmpname) // '_sfnnuc1'
	    long_name = trim(tmpname) // ' modal_aero new particle nucleation column tendency'
	    call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if ( history_aerosol ) then
               call add_default( fieldname, 1, ' ' )
            endif
	    if ( masterproc ) write(*,'(3(a,2x))') &
		'modal_aero_newnuc_init addfld', fieldname, unit
	end do ! l = ...


      return
      end subroutine modal_aero_newnuc_cam_init

  end module modal_aero_newnuc_cam
