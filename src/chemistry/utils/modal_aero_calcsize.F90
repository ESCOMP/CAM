! Portable science module for modal aerosol size calculation.
! RCE 07.04.13:  Adapted from MIRAGE2 code
module modal_aero_calcsize
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: modal_aero_calcsize_run
  public :: modal_aero_calcsize_diag_run
  public :: modal_aero_calcdry_run

  integer, public, parameter :: calcsize_nsrflx = 4

contains

subroutine modal_aero_calcsize_run( &
   ncol, pver, deltat, top_lev, &
   ntot_amode, nspec_amode, nspec_max, &
   dgnum_amode, dgnumlo_amode, dgnumhi_amode, &
   alnsg_amode, voltonumb_amode, voltonumblo_amode, voltonumbhi_amode, &
   specdens_amode, &
   mprognum_amode, &
   modeptr_aitken, modeptr_accum, &
   lmassptr_amode, numptr_amode, &
   lmassptrcw_amode, numptrcw_amode, &
   pdel, &
   gravit, pi, &
   num_q, &
   q, q_cw, &
   do_adjust, do_aitacc_transfer, &
   npair_renamexf, &
   nspecfrm_renamexf, &
   modefrm_renamexf, modetoo_renamexf, &
   lspecfrma_renamexf, lspectooa_renamexf, &
   lspecfrmc_renamexf, lspectooc_renamexf, &
   dgncur_a, &
   dqdt, dqdt_cw, &
   dotend, dotend_cw, &
   qsrflx, &
   errmsg, errflg)

   !-----------------------------------------------------------------------
   !
   ! Calculates aerosol size distribution parameters
   !    mprognum_amode >  0
   !       calculate Dgnum from mass, number, and fixed sigmag
   !    mprognum_amode <= 0
   !       calculate number from mass, fixed Dgnum, and fixed sigmag
   !
   ! Also (optionally) adjusts prognostic number to
   !    be within bounds determined by mass, Dgnum bounds, and sigma bounds
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   ! Grid and time arguments
   integer,  intent(in) :: ncol              ! number of columns
   integer,  intent(in) :: pver              ! number of vertical levels
   real(r8), intent(in) :: deltat            ! model time-step size (s)
   integer,  intent(in) :: top_lev           ! top level for aerosol calculations

   ! Mode dimension arguments
   integer,  intent(in) :: ntot_amode        ! total number of aerosol modes
   integer,  intent(in) :: nspec_amode(:)    ! number of species per mode (ntot_amode)
   integer,  intent(in) :: nspec_max         ! max number of species in any mode

   ! Mode geometry
   real(r8), intent(in) :: dgnum_amode(:)    ! default geometric mean diameter (ntot_amode)
   real(r8), intent(in) :: dgnumlo_amode(:)  ! lower bound dgnum (ntot_amode)
   real(r8), intent(in) :: dgnumhi_amode(:)  ! upper bound dgnum (ntot_amode)

   ! Derived mode quantities
   real(r8), intent(in) :: alnsg_amode(:)    ! ln(sigmag) for each mode (ntot_amode)
   real(r8), intent(in) :: voltonumb_amode(:)    ! volume-to-number (ntot_amode)
   real(r8), intent(in) :: voltonumblo_amode(:)  ! vol-to-num at dgnumlo (ntot_amode)
   real(r8), intent(in) :: voltonumbhi_amode(:)  ! vol-to-num at dgnumhi (ntot_amode)

   ! Species densities
   real(r8), intent(in) :: specdens_amode(:,:)   ! species densities (nspec_max,ntot_amode)

   ! Prognostic number flags
   integer,  intent(in) :: mprognum_amode(:)     ! prognostic number flag (ntot_amode)

   ! Mode pointers
   integer,  intent(in) :: modeptr_aitken        ! index of aitken mode
   integer,  intent(in) :: modeptr_accum          ! index of accumulation mode

   ! Species-to-q index maps (interstitial)
   integer,  intent(in) :: lmassptr_amode(:,:)   ! mass species pointer (nspec_max,ntot_amode)
   integer,  intent(in) :: numptr_amode(:)       ! number species pointer (ntot_amode)

   ! Species-to-q_cw index maps (cloud-borne)
   integer,  intent(in) :: lmassptrcw_amode(:,:) ! cloud-borne mass species pointer (nspec_max,ntot_amode)
   integer,  intent(in) :: numptrcw_amode(:)     ! cloud-borne number species pointer (ntot_amode)

   ! Atmospheric state
   real(r8), intent(in) :: pdel(:,:)         ! pressure thickness (ncol,pver)

   ! Physical constants
   real(r8), intent(in) :: gravit            ! gravitational acceleration
   real(r8), intent(in) :: pi                ! pi

   ! Species arrays
   integer,  intent(in) :: num_q             ! number of species (dimension of q/q_cw)
   real(r8), intent(in) :: q(:,:,:)          ! interstitial species (ncol,pver,num_q)
   real(r8), intent(in) :: q_cw(:,:,:)       ! cloud-borne species (ncol,pver,num_q)

   ! Control flags
   logical,  intent(in) :: do_adjust         ! adjust number to size bounds
   logical,  intent(in) :: do_aitacc_transfer ! aitken<-->accum transfer

   ! Rename transfer data (only used when do_aitacc_transfer=.true.)
   integer,  intent(in) :: npair_renamexf              ! number of rename pairs
   integer,  intent(in) :: nspecfrm_renamexf(:)        ! species count per pair
   integer,  intent(in) :: modefrm_renamexf(:)         ! from mode index per pair
   integer,  intent(in) :: modetoo_renamexf(:)         ! to mode index per pair
   integer,  intent(in) :: lspecfrma_renamexf(:,:)     ! interstitial from species indices
   integer,  intent(in) :: lspectooa_renamexf(:,:)     ! interstitial to species indices
   integer,  intent(in) :: lspecfrmc_renamexf(:,:)     ! cloud-borne from species indices
   integer,  intent(in) :: lspectooc_renamexf(:,:)     ! cloud-borne to species indices

   ! Outputs
   real(r8), intent(inout) :: dgncur_a(:,:,:)  ! dry diameter (ncol,pver,ntot_amode)
   real(r8), intent(out)   :: dqdt(:,:,:)      ! interstitial tendencies (ncol,pver,num_q)
   real(r8), intent(out)   :: dqdt_cw(:,:,:)   ! cloud-borne tendencies (ncol,pver,num_q)
   logical,  intent(out)   :: dotend(:)        ! which species have interstitial tendencies (num_q)
   logical,  intent(out)   :: dotend_cw(:)     ! which species have cloud-borne tendencies (num_q)
   real(r8), intent(out)   :: qsrflx(:,:,:,:)  ! diagnostic flux (ncol,num_q,calcsize_nsrflx,2)

   ! CCPP error reporting
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! local variables
   integer  :: i, ipair, iq
   integer  :: ixfer_acc2ait, ixfer_ait2acc
   integer  :: j, jac, jsrflx, k
   integer  :: l, l1, la, lc, lna, lnc, lsfrm, lstoo
   integer  :: n, nacc, nait

   logical  :: noxf_acc2ait(nspec_max)

   real(r8), parameter :: third = 1.0_r8/3.0_r8
   real(r8) :: delnum_a2, delnum_c2            !  work variables
   real(r8) :: delnum_a3, delnum_c3, delnum_t3 !  work variables
   real(r8) :: deltatinv                     ! 1/deltat
   real(r8) :: dgncur_c(ncol,pver,ntot_amode)
   real(r8) :: dgnyy, dgnxx                  ! dgnumlo/hi of current mode
   real(r8) :: drv_a, drv_c, drv_t           ! dry volume (cm3/mol_air)
   real(r8) :: drv_t0
   real(r8) :: drv_a_noxf, drv_c_noxf, drv_t_noxf
   real(r8) :: drv_a_acc, drv_c_acc
   real(r8) :: drv_a_accsv(ncol,pver), drv_c_accsv(ncol,pver)
   real(r8) :: drv_a_aitsv(ncol,pver), drv_c_aitsv(ncol,pver)
   real(r8) :: drv_a_sv(ncol,pver,ntot_amode), drv_c_sv(ncol,pver,ntot_amode)
   real(r8) :: dryvol_a(ncol,pver)          ! interstital aerosol dry
   ! volume (cm^3/mol_air)
   real(r8) :: dryvol_c(ncol,pver)          ! activated aerosol dry volume
   real(r8) :: duma, dumb, dumc, dumd        ! work variables
   real(r8) :: dumfac, dummwdens             ! work variables
   real(r8) :: frelaxadj                     ! relaxation factor applied
   ! to size bounds
   real(r8) :: fracadj                       ! deltat/tadj
   real(r8) :: num_a0, num_c0, num_t0        ! initial number (#/mol_air)
   real(r8) :: num_a1, num_c1                ! working number (#/mol_air)
   real(r8) :: num_a2, num_c2, num_t2        ! working number (#/mol_air)
   real(r8) :: num_a, num_c, num_t           ! final number (#/mol_air)
   real(r8) :: num_t_noxf
   real(r8) :: numbnd                        ! bounded number
   real(r8) :: num_a_acc, num_c_acc
   real(r8) :: num_a_accsv(ncol,pver), num_c_accsv(ncol,pver)
   real(r8) :: num_a_aitsv(ncol,pver), num_c_aitsv(ncol,pver)
   real(r8) :: num_a_sv(ncol,pver,ntot_amode), num_c_sv(ncol,pver,ntot_amode)
   real(r8) :: pdel_fac                      !
   real(r8) :: tadj                          ! adjustment time scale
   real(r8) :: tadjinv                       ! 1/tadj
   real(r8) :: v2ncur_a(ncol,pver,ntot_amode)
   real(r8) :: v2ncur_c(ncol,pver,ntot_amode)
   real(r8) :: v2nyy, v2nxx, v2nzz           ! voltonumblo/hi of current mode
   real(r8) :: v2nyyrl, v2nxxrl              ! relaxed voltonumblo/hi
   real(r8) :: xfercoef
   real(r8) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
   real(r8) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
   real(r8) :: xferfrac_num_acc2ait, xferfrac_vol_acc2ait
   real(r8) :: xferfrac_num_ait2acc, xferfrac_vol_ait2acc
   real(r8) :: xfertend, xfertend_num(2,2)

   integer  :: ixfer_acc2ait_sv(ncol,pver), ixfer_ait2acc_sv(ncol,pver)
   !-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   ! Initialize outputs.
   !
   ! Note on dgncur_a: it is intent(inout) and is backed by the pbuf (CAM) or
   ! non-advected constituent (SIMA) field; levels above top_lev are not modified here.
   ! It is initialized to zero or a copy of the per-mode diameter when calcsize_diag
   ! called by radiation for diagnostic lists when that is active.
   ! So it should not be initialized as part of the other (pure out) variables here.
   dotend(:) = .false.
   dotend_cw(:) = .false.
   dqdt(:,:,:) = 0.0_r8
   dqdt_cw(:,:,:) = 0.0_r8
   qsrflx(:,:,:,:) = 0.0_r8

   nait = modeptr_aitken
   nacc = modeptr_accum

   deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
   ! tadj = adjustment time scale for number, surface when they are prognosed
   !           currently set to deltat
   tadj = deltat
   tadj = 86400
   tadj = max( tadj, deltat )
   tadjinv = 1.0_r8/(tadj*(1.0_r8 + 1.0e-15_r8))
   fracadj = deltat*tadjinv
   fracadj = max( 0.0_r8, min( 1.0_r8, fracadj ) )


   !
   !
   ! the "do 40000" loop does the original (pre jan-2006)
   !   number adjustment, one mode at a time
   ! this artificially adjusts number when mean particle size is too large
   !   or too small
   !
   !
   do n = 1, ntot_amode
      ! initialize all parameters to the default values for the mode
      do k=top_lev,pver
         do i=1,ncol
            !    sgcur_a(i,k,n) = sigmag_amode(n)
            !    sgcur_c(i,k,n) = sigmag_amode(n)
            dgncur_a(i,k,n) = dgnum_amode(n)
            dgncur_c(i,k,n) = dgnum_amode(n)
            v2ncur_a(i,k,n) = voltonumb_amode(n)
            v2ncur_c(i,k,n) = voltonumb_amode(n)
            dryvol_a(i,k) = 0.0_r8
            dryvol_c(i,k) = 0.0_r8
         end do
      end do

      ! compute dry volume mixrats =
      !      sum_over_components{ component_mass mixrat / density }
      do l1 = 1, nspec_amode(n)
         ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
         dummwdens = 1.0_r8 / specdens_amode(l1,n)
         la = lmassptr_amode(l1,n)
         do k=top_lev,pver
            do i=1,ncol
               dryvol_a(i,k) = dryvol_a(i,k)    &
                  + max(0.0_r8,q(i,k,la))*dummwdens
            end do
         end do

         lc = lmassptrcw_amode(l1,n)
         do k=top_lev,pver
            do i=1,ncol
               dryvol_c(i,k) = dryvol_c(i,k)    &
                  + max(0.0_r8,q_cw(i,k,lc))*dummwdens
            end do
         end do
      end do

      ! set "short-hand" number pointers
      lna = numptr_amode(n)
      lnc = numptrcw_amode(n)


      ! go to section for appropriate number/surface diagnosed/prognosed options
      if (mprognum_amode(n) <= 0) then

         ! option 1 -- number diagnosed (fixed dgnum and sigmag)
         !    compute number tendencies that will bring numbers to their
         !    current diagnosed values
         !
         if (lna > 0) then
            dotend(lna) = .true.
            do k=top_lev,pver
               do i=1,ncol
                  dqdt(i,k,lna) = (dryvol_a(i,k)*voltonumb_amode(n)   &
                     - q(i,k,lna)) * deltatinv
               end do
            end do
         end if
         if (lnc > 0) then
            dotend_cw(lnc) = .true.
            do k=top_lev,pver
               do i=1,ncol
                  dqdt_cw(i,k,lnc) = (dryvol_c(i,k)*voltonumb_amode(n)   &
                     - q_cw(i,k,lnc)) * deltatinv
               end do
            end do
         end if
      else
         !
         ! option 2 -- number prognosed (variable dgnum, fixed sigmag)
         !       Compute number tendencies to adjust numbers if they are outside
         !    the limits determined by current volume and dgnumlo/hi
         !       The interstitial and activated aerosol fractions can, at times,
         !    be the lower or upper tail of the "total" distribution.  Thus they
         !    can be expected to have a greater range of size parameters than
         !    what is specified for the total distribution (via dgnumlo/hi)
         !       When both the interstitial and activated dry volumes are positive,
         !    the adjustment strategy is to (1) adjust the interstitial and activated
         !    numbers towards relaxed bounds, then (2) adjust the total/combined
         !    number towards the primary bounds.
         !
         ! note
         !    v2nyy = voltonumblo_amode is proportional to dgnumlo**(-3),
         !            and produces the maximum allowed number for a given volume
         !    v2nxx = voltonumbhi_amode is proportional to dgnumhi**(-3),
         !            and produces the minimum allowed number for a given volume
         !    v2nxxrl and v2nyyrl are their "relaxed" equivalents.
         !            Setting frelaxadj=27=3**3 means that
         !            dgnumlo_relaxed = dgnumlo/3 and dgnumhi_relaxed = dgnumhi*3
         !
         ! if do_aitacc_transfer is .true., then
         !     for n=nacc, multiply v2nyy by 1.0e6 to effectively turn off the
         !         adjustment when number is too big (size is too small)
         !     for n=nait, divide   v2nxx by 1.0e6 to effectively turn off the
         !         adjustment when number is too small (size is too big)
         !OLD  however, do not change the v2nyyrl/v2nxxrl so that
         !OLD      the interstitial<-->activated adjustment is not changed
         !NEW  also change the v2nyyrl/v2nxxrl so that
         !NEW      the interstitial<-->activated adjustment is turned off
         !
      end if
      frelaxadj = 27.0_r8
      dumfac = exp(4.5_r8*alnsg_amode(n)**2)*pi/6.0_r8
      v2nxx = voltonumbhi_amode(n)
      v2nyy = voltonumblo_amode(n)
      v2nxxrl = v2nxx/frelaxadj
      v2nyyrl = v2nyy*frelaxadj
      dgnxx = dgnumhi_amode(n)
      dgnyy = dgnumlo_amode(n)
      if ( do_aitacc_transfer ) then
         if (n == nait) v2nxx = v2nxx/1.0e6_r8
         if (n == nacc) v2nyy = v2nyy*1.0e6_r8
         v2nxxrl = v2nxx/frelaxadj   ! NEW
         v2nyyrl = v2nyy*frelaxadj   ! NEW
      end if

      if (do_adjust) then
         dotend(lna) = .true.
         dotend_cw(lnc) = .true.
      end if

      do  k = top_lev, pver
         do  i = 1, ncol

            drv_a = dryvol_a(i,k)
            num_a0 = q(i,k,lna)
            num_a = max( 0.0_r8, num_a0 )
            drv_c = dryvol_c(i,k)
            num_c0 = q_cw(i,k,lnc)
            num_c = max( 0.0_r8, num_c0 )

            if (do_adjust) then

               !
               ! do number adjustment for interstitial and activated particles
               !    adjustments that (1) make numbers non-negative or (2) make numbers
               !       zero when volume is zero are applied over time-scale deltat
               !    adjustments that bring numbers to within specified bounds are
               !       applied over time-scale tadj
               !
               if ((drv_a <= 0.0_r8) .and. (drv_c <= 0.0_r8)) then
                  ! both interstitial and activated volumes are zero
                  ! adjust both numbers to zero
                  num_a = 0.0_r8
                  dqdt(i,k,lna) = -num_a0*deltatinv
                  num_c = 0.0_r8
                  dqdt_cw(i,k,lnc) = -num_c0*deltatinv
               else if (drv_c <= 0.0_r8) then
                  ! activated volume is zero, so interstitial number/volume == total/combined
                  ! apply step 1 and 3, but skip the relaxed adjustment (step 2, see below)
                  num_c = 0.0_r8
                  dqdt_cw(i,k,lnc) = -num_c0*deltatinv
                  num_a1 = num_a
                  numbnd = max( drv_a*v2nxx, min( drv_a*v2nyy, num_a1 ) )
                  num_a  = num_a1 + (numbnd - num_a1)*fracadj
                  dqdt(i,k,lna) = (num_a - num_a0)*deltatinv

               else if (drv_a <= 0.0_r8) then
                  ! interstitial volume is zero, treat similar to above
                  num_a = 0.0_r8
                  dqdt(i,k,lna) = -num_a0*deltatinv
                  num_c1 = num_c
                  numbnd = max( drv_c*v2nxx, min( drv_c*v2nyy, num_c1 ) )
                  num_c  = num_c1 + (numbnd - num_c1)*fracadj
                  dqdt_cw(i,k,lnc) = (num_c - num_c0)*deltatinv
               else
                  ! both volumes are positive
                  ! apply 3 adjustment steps
                  ! step1:  num_a,c0 --> num_a,c1 forces non-negative values
                  num_a1 = num_a
                  num_c1 = num_c
                  ! step2:  num_a,c1 --> num_a,c2 applies relaxed bounds to the interstitial
                  !    and activated number (individually)
                  !    if only only a or c changes, adjust the other in the opposite direction
                  !    as much as possible to conserve a+c
                  numbnd = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl, num_a1 ) )
                  delnum_a2 = (numbnd - num_a1)*fracadj
                  num_a2 = num_a1 + delnum_a2
                  numbnd = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl, num_c1 ) )
                  delnum_c2 = (numbnd - num_c1)*fracadj
                  num_c2 = num_c1 + delnum_c2
                  if ((delnum_a2 == 0.0_r8) .and. (delnum_c2 /= 0.0_r8)) then
                     num_a2 = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl,   &
                        num_a1-delnum_c2 ) )
                  else if ((delnum_a2 /= 0.0_r8) .and. (delnum_c2 == 0.0_r8)) then
                     num_c2 = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl,   &
                        num_c1-delnum_a2 ) )
                  end if
                  ! step3:  num_a,c2 --> num_a,c3 applies stricter bounds to the
                  !    combined/total number
                  drv_t = drv_a + drv_c
                  num_t2 = num_a2 + num_c2
                  delnum_a3 = 0.0_r8
                  delnum_c3 = 0.0_r8
                  if (num_t2 < drv_t*v2nxx) then
                     delnum_t3 = (drv_t*v2nxx - num_t2)*fracadj
                     ! if you are here then (num_a2 < drv_a*v2nxx) and/or
                     !                      (num_c2 < drv_c*v2nxx) must be true
                     if ((num_a2 < drv_a*v2nxx) .and. (num_c2 < drv_c*v2nxx)) then
                        delnum_a3 = delnum_t3*(num_a2/num_t2)
                        delnum_c3 = delnum_t3*(num_c2/num_t2)
                     else if (num_c2 < drv_c*v2nxx) then
                        delnum_c3 = delnum_t3
                     else if (num_a2 < drv_a*v2nxx) then
                        delnum_a3 = delnum_t3
                     end if
                  else if (num_t2 > drv_t*v2nyy) then
                     delnum_t3 = (drv_t*v2nyy - num_t2)*fracadj
                     ! if you are here then (num_a2 > drv_a*v2nyy) and/or
                     !                      (num_c2 > drv_c*v2nyy) must be true
                     if ((num_a2 > drv_a*v2nyy) .and. (num_c2 > drv_c*v2nyy)) then
                        delnum_a3 = delnum_t3*(num_a2/num_t2)
                        delnum_c3 = delnum_t3*(num_c2/num_t2)
                     else if (num_c2 > drv_c*v2nyy) then
                        delnum_c3 = delnum_t3
                     else if (num_a2 > drv_a*v2nyy) then
                        delnum_a3 = delnum_t3
                     end if
                  end if
                  num_a = num_a2 + delnum_a3
                  dqdt(i,k,lna) = (num_a - num_a0)*deltatinv
                  num_c = num_c2 + delnum_c3
                  dqdt_cw(i,k,lnc) = (num_c - num_c0)*deltatinv
               end if

            end if ! do_adjust

            !
            ! now compute current dgn and v2n
            !
            if (drv_a > 0.0_r8) then
               if (num_a <= drv_a*v2nxx) then
                  dgncur_a(i,k,n) = dgnxx
                  v2ncur_a(i,k,n) = v2nxx
               else if (num_a >= drv_a*v2nyy) then
                  dgncur_a(i,k,n) = dgnyy
                  v2ncur_a(i,k,n) = v2nyy
               else
                  dgncur_a(i,k,n) = (drv_a/(dumfac*num_a))**third
                  v2ncur_a(i,k,n) = num_a/drv_a
               end if
            end if
            pdel_fac = pdel(i,k)/gravit   ! = rho*dz
            jac = 1
            qsrflx(i,lna,1,jac) = qsrflx(i,lna,1,jac) + max(0.0_r8,dqdt(i,k,lna))*pdel_fac
            qsrflx(i,lna,2,jac) = qsrflx(i,lna,2,jac) + min(0.0_r8,dqdt(i,k,lna))*pdel_fac

            if (drv_c > 0.0_r8) then
               if (num_c <= drv_c*v2nxx) then
                  dgncur_c(i,k,n) = dgnumhi_amode(n)
                  v2ncur_c(i,k,n) = v2nxx
               else if (num_c >= drv_c*v2nyy) then
                  dgncur_c(i,k,n) = dgnumlo_amode(n)
                  v2ncur_c(i,k,n) = v2nyy
               else
                  dgncur_c(i,k,n) = (drv_c/(dumfac*num_c))**third
                  v2ncur_c(i,k,n) = num_c/drv_c
               end if
            end if
            jac = 2
            qsrflx(i,lnc,1,jac) = qsrflx(i,lnc,1,jac) + max(0.0_r8,dqdt_cw(i,k,lnc))*pdel_fac
            qsrflx(i,lnc,2,jac) = qsrflx(i,lnc,2,jac) + min(0.0_r8,dqdt_cw(i,k,lnc))*pdel_fac


            ! save number and dryvol for aitken <--> accum renaming
            if ( do_aitacc_transfer ) then
               if (n == nait) then
                  drv_a_aitsv(i,k) = drv_a
                  num_a_aitsv(i,k) = num_a
                  drv_c_aitsv(i,k) = drv_c
                  num_c_aitsv(i,k) = num_c
               else if (n == nacc) then
                  drv_a_accsv(i,k) = drv_a
                  num_a_accsv(i,k) = num_a
                  drv_c_accsv(i,k) = drv_c
                  num_c_accsv(i,k) = num_c
               end if
            end if
            drv_a_sv(i,k,n) = drv_a
            num_a_sv(i,k,n) = num_a
            drv_c_sv(i,k,n) = drv_c
            num_c_sv(i,k,n) = num_c

         end do
      end do


      !
      ! option 3 -- number and surface prognosed (variable dgnum and sigmag)
      !             this is not implemented
      !
   end do  ! do n = 1, ntot_amode


   !
   !
   ! the following section
   !    does aitken <--> accum mode transfer
   !
   ! when the aitken mode mean size is too big, the largest
   !    aitken particles are transferred into the accum mode
   !    to reduce the aitken mode mean size
   ! when the accum mode mean size is too small, the smallest
   !    accum particles are transferred into the aitken mode
   !    to increase the accum mode mean size
   !
   !
   ixfer_ait2acc_sv(:,:) = 0
   ixfer_acc2ait_sv(:,:) = 0
   if ( do_aitacc_transfer ) then

      if (npair_renamexf .le. 0) then
         errmsg = 'modal_aero_calcsize_run error -- npair_renamexf <= 0'
         errflg = 1
         return
      end if

      ! check that renaming ipair=1 is aitken-->accum
      ipair = 1
      if ((modefrm_renamexf(ipair) .ne. nait) .or.   &
         (modetoo_renamexf(ipair) .ne. nacc)) then
         errmsg = 'modal_aero_calcsize_run error -- modefrm/too_renamexf(1) are wrong'
         errflg = 1
         return
      end if

      ! set dotend() for species that will be transferred
      do iq = 1, nspecfrm_renamexf(ipair)
         lsfrm = lspecfrma_renamexf(iq,ipair)
         lstoo = lspectooa_renamexf(iq,ipair)
         if ((lsfrm > 0) .and. (lstoo > 0)) then
            dotend(lsfrm) = .true.
            dotend(lstoo) = .true.
         end if
         lsfrm = lspecfrmc_renamexf(iq,ipair)
         lstoo = lspectooc_renamexf(iq,ipair)
         if ((lsfrm > 0) .and. (lstoo > 0)) then
            dotend_cw(lsfrm) = .true.
            dotend_cw(lstoo) = .true.
         end if
      end do

      ! identify accum species cannot be transferred to aitken mode
      noxf_acc2ait(:) = .true.
      do l1 = 1, nspec_amode(nacc)
         la = lmassptr_amode(l1,nacc)
         do iq = 1, nspecfrm_renamexf(ipair)
            if (lspectooa_renamexf(iq,ipair) == la) then
               noxf_acc2ait(l1) = .false.
            end if
         end do
      end do

      ! v2nzz is voltonumb at the "geometrically-defined" mid-point
      ! between the aitken and accum modes
      v2nzz = sqrt(voltonumb_amode(nait)*voltonumb_amode(nacc))

      ! loop over columns and levels
      do  k = top_lev, pver
         do  i = 1, ncol

            pdel_fac = pdel(i,k)/gravit   ! = rho*dz
            xfertend_num(:,:) = 0.0_r8

            ! compute aitken --> accum transfer rates
            ixfer_ait2acc = 0
            xfercoef_num_ait2acc = 0.0_r8
            xfercoef_vol_ait2acc = 0.0_r8

            drv_t = drv_a_aitsv(i,k) + drv_c_aitsv(i,k)
            num_t = num_a_aitsv(i,k) + num_c_aitsv(i,k)
            if (drv_t > 0.0_r8) then
               if (num_t < drv_t*v2nzz) then
                  ixfer_ait2acc = 1
                  if (num_t < drv_t*voltonumb_amode(nacc)) then
                     xferfrac_num_ait2acc = 1.0_r8
                     xferfrac_vol_ait2acc = 1.0_r8
                  else
                     xferfrac_vol_ait2acc = ((num_t/drv_t) - v2nzz)/   &
                        (voltonumb_amode(nacc) - v2nzz)
                     xferfrac_num_ait2acc = xferfrac_vol_ait2acc*   &
                        (drv_t*voltonumb_amode(nacc)/num_t)
                     if ((xferfrac_num_ait2acc <= 0.0_r8) .or.   &
                        (xferfrac_vol_ait2acc <= 0.0_r8)) then
                        xferfrac_num_ait2acc = 0.0_r8
                        xferfrac_vol_ait2acc = 0.0_r8
                     else if ((xferfrac_num_ait2acc >= 1.0_r8) .or.   &
                        (xferfrac_vol_ait2acc >= 1.0_r8)) then
                        xferfrac_num_ait2acc = 1.0_r8
                        xferfrac_vol_ait2acc = 1.0_r8
                     end if
                  end if
                  xfercoef_num_ait2acc = xferfrac_num_ait2acc*tadjinv
                  xfercoef_vol_ait2acc = xferfrac_vol_ait2acc*tadjinv
                  xfertend_num(1,1) = num_a_aitsv(i,k)*xfercoef_num_ait2acc
                  xfertend_num(1,2) = num_c_aitsv(i,k)*xfercoef_num_ait2acc
               end if
            end if

            ! compute accum --> aitken transfer rates
            ! accum may have some species (seasalt, dust, poa, lll) that are
            !    not in aitken mode
            ! so first divide the accum drv & num into not-transferred (noxf) species
            !    and transferred species, and use the transferred-species
            !    portion in what follows
            ixfer_acc2ait = 0
            xfercoef_num_acc2ait = 0.0_r8
            xfercoef_vol_acc2ait = 0.0_r8

            drv_t = drv_a_accsv(i,k) + drv_c_accsv(i,k)
            num_t = num_a_accsv(i,k) + num_c_accsv(i,k)
            drv_a_noxf = 0.0_r8
            drv_c_noxf = 0.0_r8
            if (drv_t > 0.0_r8) then
               if (num_t > drv_t*v2nzz) then
                  do l1 = 1, nspec_amode(nacc)

                     if ( noxf_acc2ait(l1) ) then
                        ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
                        dummwdens = 1.0_r8 / specdens_amode(l1,nacc)
                        la = lmassptr_amode(l1,nacc)
                        drv_a_noxf = drv_a_noxf    &
                           + max(0.0_r8,q(i,k,la))*dummwdens
                        lc = lmassptrcw_amode(l1,nacc)

                        drv_c_noxf = drv_c_noxf    &
                           + max(0.0_r8,q_cw(i,k,lc))*dummwdens
                     end if
                  end do
                  drv_t_noxf = drv_a_noxf + drv_c_noxf
                  num_t_noxf = drv_t_noxf*voltonumblo_amode(nacc)
                  num_t0 = num_t
                  drv_t0 = drv_t
                  num_t = max( 0.0_r8, num_t - num_t_noxf )
                  drv_t = max( 0.0_r8, drv_t - drv_t_noxf )
               end if
            end if

            if (drv_t > 0.0_r8) then
               if (num_t > drv_t*v2nzz) then
                  ixfer_acc2ait = 1
                  if (num_t > drv_t*voltonumb_amode(nait)) then
                     xferfrac_num_acc2ait = 1.0_r8
                     xferfrac_vol_acc2ait = 1.0_r8
                  else
                     xferfrac_vol_acc2ait = ((num_t/drv_t) - v2nzz)/   &
                        (voltonumb_amode(nait) - v2nzz)
                     xferfrac_num_acc2ait = xferfrac_vol_acc2ait*   &
                        (drv_t*voltonumb_amode(nait)/num_t)
                     if ((xferfrac_num_acc2ait <= 0.0_r8) .or.   &
                        (xferfrac_vol_acc2ait <= 0.0_r8)) then
                        xferfrac_num_acc2ait = 0.0_r8
                        xferfrac_vol_acc2ait = 0.0_r8
                     else if ((xferfrac_num_acc2ait >= 1.0_r8) .or.   &
                        (xferfrac_vol_acc2ait >= 1.0_r8)) then
                        xferfrac_num_acc2ait = 1.0_r8
                        xferfrac_vol_acc2ait = 1.0_r8
                     end if
                  end if
                  duma = 1.0e-37_r8
                  xferfrac_num_acc2ait = xferfrac_num_acc2ait*   &
                     num_t/max( duma, num_t0 )
                  xfercoef_num_acc2ait = xferfrac_num_acc2ait*tadjinv
                  xfercoef_vol_acc2ait = xferfrac_vol_acc2ait*tadjinv
                  xfertend_num(2,1) = num_a_accsv(i,k)*xfercoef_num_acc2ait
                  xfertend_num(2,2) = num_c_accsv(i,k)*xfercoef_num_acc2ait
               end if
            end if

            ! jump to end-of-loop if no transfer is needed at current i,k
            if (ixfer_ait2acc+ixfer_acc2ait > 0) then
               ixfer_ait2acc_sv(i,k) = ixfer_ait2acc
               ixfer_acc2ait_sv(i,k) = ixfer_acc2ait

               !
               ! compute new dgncur & v2ncur for aitken & accum modes
               !
               ! currently inactive
               do n = nait, nacc, (nacc-nait)
                  if (n .eq. nait) then
                     duma = (xfertend_num(1,1) - xfertend_num(2,1))*deltat
                     num_a     = max( 0.0_r8, num_a_aitsv(i,k) - duma )
                     num_a_acc = max( 0.0_r8, num_a_accsv(i,k) + duma )
                     duma = (drv_a_aitsv(i,k)*xfercoef_vol_ait2acc -   &
                        (drv_a_accsv(i,k)-drv_a_noxf)*xfercoef_vol_acc2ait)*deltat
                     drv_a     = max( 0.0_r8, drv_a_aitsv(i,k) - duma )
                     drv_a_acc = max( 0.0_r8, drv_a_accsv(i,k) + duma )
                     duma = (xfertend_num(1,2) - xfertend_num(2,2))*deltat
                     num_c     = max( 0.0_r8, num_c_aitsv(i,k) - duma )
                     num_c_acc = max( 0.0_r8, num_c_accsv(i,k) + duma )
                     duma = (drv_c_aitsv(i,k)*xfercoef_vol_ait2acc -   &
                        (drv_c_accsv(i,k)-drv_c_noxf)*xfercoef_vol_acc2ait)*deltat
                     drv_c     = max( 0.0_r8, drv_c_aitsv(i,k) - duma )
                     drv_c_acc = max( 0.0_r8, drv_c_accsv(i,k) + duma )
                  else
                     num_a = num_a_acc
                     drv_a = drv_a_acc
                     num_c = num_c_acc
                     drv_c = drv_c_acc
                  end if

                  if (drv_a > 0.0_r8) then
                     if (num_a <= drv_a*voltonumbhi_amode(n)) then
                        dgncur_a(i,k,n) = dgnumhi_amode(n)
                        v2ncur_a(i,k,n) = voltonumbhi_amode(n)
                     else if (num_a >= drv_a*voltonumblo_amode(n)) then
                        dgncur_a(i,k,n) = dgnumlo_amode(n)
                        v2ncur_a(i,k,n) = voltonumblo_amode(n)
                     else
                        dgncur_a(i,k,n) = (drv_a/(dumfac*num_a))**third
                        v2ncur_a(i,k,n) = num_a/drv_a
                     end if
                  else
                     dgncur_a(i,k,n) = dgnum_amode(n)
                     v2ncur_a(i,k,n) = voltonumb_amode(n)
                  end if

                  if (drv_c > 0.0_r8) then
                     if (num_c <= drv_c*voltonumbhi_amode(n)) then
                        dgncur_c(i,k,n) = dgnumhi_amode(n)
                        v2ncur_c(i,k,n) = voltonumbhi_amode(n)
                     else if (num_c >= drv_c*voltonumblo_amode(n)) then
                        dgncur_c(i,k,n) = dgnumlo_amode(n)
                        v2ncur_c(i,k,n) = voltonumblo_amode(n)
                     else
                        dgncur_c(i,k,n) = (drv_c/(dumfac*num_c))**third
                        v2ncur_c(i,k,n) = num_c/drv_c
                     end if
                  else
                     dgncur_c(i,k,n) = dgnum_amode(n)
                     v2ncur_c(i,k,n) = voltonumb_amode(n)
                  end if

               end do


               !
               ! compute tendency amounts for aitken <--> accum transfer
               !

               ! j=1 does aitken-->accum; j=2 does accum-->aitken
               do  j = 1, 2

                  if ((j .eq. 1 .and. ixfer_ait2acc > 0) .or. &
                     (j .eq. 2 .and. ixfer_acc2ait > 0)) then

                     jsrflx = j+2
                     if (j .eq. 1) then
                        xfercoef = xfercoef_vol_ait2acc
                     else
                        xfercoef = xfercoef_vol_acc2ait
                     end if

                     do  iq = 1, nspecfrm_renamexf(ipair)

                        ! jac=1 does interstitial ("_a"); jac=2 does activated ("_c");
                        do  jac = 1, 2

                           ! the lspecfrma_renamexf (and lspecfrmc_renamexf) are aitken species
                           ! the lspectooa_renamexf (and lspectooc_renamexf) are accum  species
                           ! for j=1, want lsfrm=aitken species, lstoo=accum  species
                           ! for j=2, want lsfrm=accum  species,  lstoo=aitken species
                           if (j .eq. 1) then
                              if (jac .eq. 1) then
                                 lsfrm = lspecfrma_renamexf(iq,ipair)
                                 lstoo = lspectooa_renamexf(iq,ipair)
                              else
                                 lsfrm = lspecfrmc_renamexf(iq,ipair)
                                 lstoo = lspectooc_renamexf(iq,ipair)
                              end if
                           else
                              if (jac .eq. 1) then
                                 lsfrm = lspectooa_renamexf(iq,ipair)
                                 lstoo = lspecfrma_renamexf(iq,ipair)
                              else
                                 lsfrm = lspectooc_renamexf(iq,ipair)
                                 lstoo = lspecfrmc_renamexf(iq,ipair)
                              end if
                           end if

                           if ((lsfrm > 0) .and. (lstoo > 0)) then
                              if (jac .eq. 1) then
                                 if (iq .eq. 1) then
                                    xfertend = xfertend_num(j,jac)
                                 else
                                    xfertend = max(0.0_r8,q(i,k,lsfrm))*xfercoef
                                 end if
                                 dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xfertend
                                 dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xfertend
                              else
                                 if (iq .eq. 1) then
                                    xfertend = xfertend_num(j,jac)
                                 else
                                    xfertend = max(0.0_r8,q_cw(i,k,lsfrm))*xfercoef
                                 end if
                                 dqdt_cw(i,k,lsfrm) = dqdt_cw(i,k,lsfrm) - xfertend
                                 dqdt_cw(i,k,lstoo) = dqdt_cw(i,k,lstoo) + xfertend
                              end if
                              qsrflx(i,lsfrm,jsrflx,jac) = qsrflx(i,lsfrm,jsrflx,jac) - xfertend*pdel_fac
                              qsrflx(i,lstoo,jsrflx,jac) = qsrflx(i,lstoo,jsrflx,jac) + xfertend*pdel_fac
                           end if

                        end do
                     end do
                  end if
               end do

            end if
         end do
      end do


   end if  !  do_aitacc_transfer
   lsfrm = -123456789   ! executable statement for debugging

end subroutine modal_aero_calcsize_run

!===============================================================================

subroutine modal_aero_calcdry_run( &
   aero_props, aero_state, &
   ncol, pver, top_lev, &
   do_strat_sulfate, &
   pi, &
   dgncur_a, &
   hygro, dryvol, dryrad, drymass, so4dryvol, naer, &
   errmsg, errflg)

!-----------------------------------------------------------------------
!
! Compute derived dry aerosol properties from mixing ratios and
! adjusted number mode diameter. Called after calcsize_run.
!
!-----------------------------------------------------------------------

   use aerosol_properties_mod, only: aerosol_properties
   use aerosol_state_mod,      only: aerosol_state

   ! Arguments
   class(aerosol_properties), intent(in) :: aero_props
   class(aerosol_state),      intent(in) :: aero_state
   integer,          intent(in)  :: ncol                    ! number of columns
   integer,          intent(in)  :: pver                    ! number of vertical levels
   integer,          intent(in)  :: top_lev                 ! top level for aerosol calculations
   logical,          intent(in)  :: do_strat_sulfate        ! use stratospheric sulfate treatment
   real(r8),         intent(in)  :: pi                      ! pi
   real(r8),         intent(in)  :: dgncur_a(:,:,:)         ! dry number mode diameter (m)

   real(r8),         intent(out) :: hygro(:,:,:)            ! volume-weighted mean hygroscopicity (--)
   real(r8),         intent(out) :: dryvol(:,:,:)           ! single-particle-mean dry volume (m3)
   real(r8),         intent(out) :: dryrad(:,:,:)           ! dry volume mean radius of aerosol (m)
   real(r8),         intent(out) :: drymass(:,:,:)          ! single-particle-mean dry mass (kg)
   real(r8),         intent(out) :: so4dryvol(:,:,:)        ! single-particle-mean so4 dry volume (m3)
   real(r8),         intent(out) :: naer(:,:,:)             ! aerosol number MR (#/kg-air)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! local variables
   real(r8), parameter :: third = 1._r8/3._r8
   real(r8) :: pi43

   integer  :: i, k, l, m
   integer  :: nmodes, nspec

   real(r8) :: specdens
   real(r8) :: spechygro, spechygro_1
   real(r8) :: sigmag
   real(r8) :: duma, dumb
   real(r8) :: alnsg

   real(r8) :: v2ncur_a
   real(r8) :: drydens               ! dry particle density  (kg/m^3)

   real(r8) :: maer(ncol, pver)
   real(r8) :: dryvolmr(ncol, pver)
   real(r8) :: so4dryvolmr(ncol, pver)

   character(len=32) :: spectype

   real(r8), pointer :: raer(:,:)
   !-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   pi43 = pi*4._r8/3._r8

   nmodes = aero_props%nbins()

   hygro(:,:,:)     = 0._r8
   so4dryvol(:,:,:) = 0._r8

   do m = 1, nmodes

      maer(:,:)        = 0._r8
      dryvolmr(:,:)    = 0._r8
      so4dryvolmr(:,:) = 0._r8

      ! get mode properties
      sigmag = exp(aero_props%alogsig(m))

      ! get mode info
      nspec = aero_props%nspecies(m)

      do l = 1, nspec

         ! get species interstitial mixing ratio ('a')
         call aero_state%get_ambient_mmr(species_ndx=l, bin_ndx=m, mmr=raer)
         call aero_props%get(m, l, density=specdens, &
                                     hygro=spechygro, spectype=spectype)

         if (l == 1) then
            ! save off these values to be used as defaults
            spechygro_1    = spechygro
         end if

         do k = top_lev, pver
            do i = 1, ncol
               duma          = raer(i,k)     ! kg/kg air
               maer(i,k)     = maer(i,k) + duma
               dumb          = duma/specdens ! m3/kg air
               dryvolmr(i,k) = dryvolmr(i,k) + dumb
               if (do_strat_sulfate .and. (trim(spectype).eq.'sulfate')) then
                  so4dryvolmr(i,k) = so4dryvolmr(i,k) + dumb
               end if
               hygro(i,k,m)  = hygro(i,k,m) + dumb*spechygro
            end do
         end do
      end do

      alnsg = log(sigmag)

      do k = top_lev, pver
         do i = 1, ncol

            if (dryvolmr(i,k) > 1.0e-30_r8) then
               hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k)
            else
               hygro(i,k,m) = spechygro_1
            end if

            ! dry aerosol properties

            v2ncur_a = 1._r8 / ( (pi/6._r8)*(dgncur_a(i,k,m)**3._r8)*exp(4.5_r8*alnsg**2._r8) )
            ! naer = aerosol number (#/kg)
            naer(i,k,m) = dryvolmr(i,k)*v2ncur_a

            ! compute mean (1 particle) dry volume and mass for each mode
            if (maer(i,k) .gt. 1.0e-31_r8) then
               drydens = maer(i,k)/dryvolmr(i,k)        ! kg/m3 aerosol
            else
               drydens = 1.0_r8
            end if
            dryvol(i,k,m)   = 1.0_r8/v2ncur_a             ! m3/particle
            drymass(i,k,m)  = drydens*dryvol(i,k,m)       ! kg/particle
            dryrad(i,k,m)   = (dryvol(i,k,m)/pi43)**third ! m
         end do    ! i = 1, ncol
      end do    ! k = top_lev, pver


      if (do_strat_sulfate) then
         do k = top_lev, pver
            do i = 1, ncol
               if (so4dryvolmr(i,k) .gt. 1.0e-31_r8) then
                  so4dryvol(i,k,m) = dryvol(i,k,m)*so4dryvolmr(i,k)/dryvolmr(i,k)
               else
                  so4dryvol(i,k,m) = 0.0_r8
               end if

            end do    ! i = 1, ncol
         end do    ! k = top_lev, pver

      end if

   end do    ! m = 1, nmodes

end subroutine modal_aero_calcdry_run

subroutine modal_aero_calcsize_diag_run( &
   aero_props, aero_state, &
   ncol, pver, top_lev, &
   pi, &
   dgncur_a, &
   errmsg, errflg)

   !-----------------------------------------------------------------------
   !
   ! Calculate aerosol size distribution parameters for a diagnostic
   ! radiation list, using only the abstract aerosol interfaces.
   ! Number is diagnosed from mass, Dgnum bounds, and fixed sigmag
   ! (the mprognum <= 0 branch of the prognostic calculation).
   !
   !-----------------------------------------------------------------------

   use aerosol_properties_mod, only: aerosol_properties
   use aerosol_state_mod,      only: aerosol_state

   ! Arguments
   class(aerosol_properties), intent(in) :: aero_props
   class(aerosol_state),      intent(in) :: aero_state
   integer,          intent(in)  :: ncol              ! number of columns
   integer,          intent(in)  :: pver              ! number of vertical levels
   integer,          intent(in)  :: top_lev           ! top level for aerosol calculations
   real(r8),         intent(in)  :: pi                ! pi
   real(r8),         intent(out) :: dgncur_a(:,:,:)   ! dry number mode diameter (m)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! local
   integer  :: i, k, l1, n
   integer  :: nmodes
   integer  :: nspec

   real(r8), parameter :: third = 1.0_r8/3.0_r8

   real(r8), pointer :: mode_num(:,:) ! mode number mixing ratio
   real(r8), pointer :: specmmr(:,:)  ! specie mmr
   real(r8)          :: specdens      ! specie density

   real(r8) :: dryvol_a(ncol,pver)    ! interstital aerosol dry volume (cm^3/mol_air)

   real(r8) :: dgnum, dgnumhi, dgnumlo
   real(r8) :: dgnyy, dgnxx           ! dgnumlo/hi of current mode
   real(r8) :: drv_a                  ! dry volume (cm3/mol_air)
   real(r8) :: dumfac, dummwdens      ! work variables
   real(r8) :: num_a0                 ! initial number (#/mol_air)
   real(r8) :: num_a                  ! final number (#/mol_air)
   real(r8) :: voltonumbhi, voltonumblo
   real(r8) :: v2nyy, v2nxx           ! voltonumblo/hi of current mode
   real(r8) :: sigmag, alnsg
   !-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   nmodes = aero_props%nbins()

   do n = 1, nmodes

      ! get mode properties
      dgnum = aero_props%dgnum(n)
      dgnumhi = aero_props%dgnumhi(n)
      dgnumlo = aero_props%dgnumlo(n)
      sigmag = exp(aero_props%alogsig(n))

      ! get mode number mixing ratio
      call aero_state%get_ambient_num(n, mode_num)

      dgncur_a(:,:,n) = dgnum
      dryvol_a(:,:) = 0.0_r8

      ! compute dry volume mixrats =
      !      sum_over_components{ component_mass mixrat / density }
      nspec = aero_props%nspecies(n)
      do l1 = 1, nspec

         call aero_state%get_ambient_mmr(species_ndx=l1, bin_ndx=n, mmr=specmmr)
         call aero_props%get(n, l1, density=specdens)

         ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
         dummwdens = 1.0_r8 / specdens

         do k=top_lev,pver
            do i=1,ncol
               dryvol_a(i,k) = dryvol_a(i,k)    &
                  + max(0.0_r8, specmmr(i,k))*dummwdens
            end do
         end do
      end do

      alnsg  = log( sigmag )
      dumfac = exp(4.5_r8*alnsg**2)*pi/6.0_r8
      voltonumblo = 1._r8 / ( (pi/6._r8)*(dgnumlo**3)*exp(4.5_r8*alnsg**2) )
      voltonumbhi = 1._r8 / ( (pi/6._r8)*(dgnumhi**3)*exp(4.5_r8*alnsg**2) )
      v2nxx = voltonumbhi
      v2nyy = voltonumblo
      dgnxx = dgnumhi
      dgnyy = dgnumlo

      do k = top_lev, pver
         do i = 1, ncol

            drv_a = dryvol_a(i,k)
            num_a0 = mode_num(i,k)
            num_a = max( 0.0_r8, num_a0 )

            if (drv_a > 0.0_r8) then
               if (num_a <= drv_a*v2nxx) then
                  dgncur_a(i,k,n) = dgnxx
               else if (num_a >= drv_a*v2nyy) then
                  dgncur_a(i,k,n) = dgnyy
               else
                  dgncur_a(i,k,n) = (drv_a/(dumfac*num_a))**third
               end if
            end if

         end do
      end do

   end do ! nmodes

end subroutine modal_aero_calcsize_diag_run

end module modal_aero_calcsize
