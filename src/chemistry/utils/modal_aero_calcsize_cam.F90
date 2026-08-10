module modal_aero_calcsize_cam
! CAM wrapper for modal_aero_calcsize.

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use physconst,        only: pi, rhoh2o, gravit

use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state, physics_ptend
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, &
                            pbuf_get_field, pbuf_add_field, pbuf_set_field, dtype_r8

use phys_control,     only: phys_getopts
use aerosol_properties_mod, only: aerosol_properties
use aerosol_state_mod, only: aerosol_state

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use cam_history,      only: addfld, add_default, fieldname_len, horiz_only, outfld
use constituents,     only: pcnst, cnst_name
use modal_aero_data,          only: modal_strat_sulfate

use ref_pres,         only: top_lev => clim_modal_aero_top_lev

#ifdef MODAL_AERO

use modal_aero_data, only: ntot_amode, nspec_amode, nspec_max, &
                           numptr_amode, &
                           alnsg_amode, &
                           voltonumbhi_amode, voltonumblo_amode, &
                           dgnum_amode, dgnumhi_amode, dgnumlo_amode

use modal_aero_data,  only: numptrcw_amode, mprognum_amode, qqcw_get_field, lmassptrcw_amode, &
           lmassptr_amode, modeptr_accum, modeptr_aitken, &
           specmw_amode, specdens_amode, voltonumb_amode, &
           cnst_name_cw

use modal_aero_rename_cam, only: lspectooa_renamexf, lspecfrma_renamexf, lspectooc_renamexf, lspecfrmc_renamexf, &
           modetoo_renamexf, nspecfrm_renamexf, npair_renamexf, modefrm_renamexf

#endif

use modal_aero_calcsize, only: modal_aero_calcsize_run, modal_aero_calcdry_run, calcsize_nsrflx

implicit none
private
save

public modal_aero_calcsize_init, modal_aero_calcsize_sub, modal_aero_calcsize_diag
public :: modal_aero_calcsize_reg

logical :: do_adjust_default
logical :: do_aitacc_transfer_default

integer :: dgnum_idx     = -1
integer :: hygro_idx     = -1
integer :: dryvol_idx    = -1
integer :: dryrad_idx    = -1
integer :: drymass_idx   = -1
integer :: so4dryvol_idx = -1
integer :: naer_idx      = -1
integer :: sulfeq_idx    = -1


!===============================================================================
contains
!===============================================================================

subroutine modal_aero_calcsize_reg()
  use radiative_aerosol, only: rad_aer_get_info

  integer :: nmodes

  call rad_aer_get_info(0, nmodes=nmodes)

  call pbuf_add_field('DGNUM', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnum_idx)

  call pbuf_add_field('HYGRO',     'phys_pkg', dtype_r8, (/pcols,pver,nmodes/), hygro_idx)
  call pbuf_add_field('DRYVOL',    'phys_pkg', dtype_r8, (/pcols,pver,nmodes/), dryvol_idx)
  call pbuf_add_field('DRYRAD',    'phys_pkg', dtype_r8, (/pcols,pver,nmodes/), dryrad_idx)
  call pbuf_add_field('DRYMASS',   'phys_pkg', dtype_r8, (/pcols,pver,nmodes/), drymass_idx)
  call pbuf_add_field('SO4DRYVOL', 'phys_pkg', dtype_r8, (/pcols,pver,nmodes/), so4dryvol_idx)
  call pbuf_add_field('NAER',      'phys_pkg', dtype_r8, (/pcols,pver,nmodes/), naer_idx)

end subroutine modal_aero_calcsize_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_calcsize_init(pbuf2d)
   use time_manager,  only: is_first_step

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   !    set do_adjust_default and do_aitacc_transfer_default flags
   !    create history fields for column tendencies associated with
   !       modal_aero_calcsize
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local
   integer  :: ipair, iq
   integer  :: jac
   integer  :: lsfrm, lstoo
   integer  :: n, nacc, nait
   logical  :: history_aerosol

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit
   !-----------------------------------------------------------------------

   call phys_getopts(history_aerosol_out=history_aerosol)

   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, dgnum_idx, 0.0_r8)
   endif

#ifndef MODAL_AERO
   do_adjust_default          = .false.
   do_aitacc_transfer_default = .false.
#else
   do_adjust_default = .true.

   nait = modeptr_aitken
   nacc = modeptr_accum
   do_aitacc_transfer_default = .false.
   if ((modeptr_aitken > 0) .and.   &
      (modeptr_accum  > 0) .and.   &
      (modeptr_aitken /= modeptr_accum)) then
      do_aitacc_transfer_default = .true.
      if (mprognum_amode(nait) <= 0) do_aitacc_transfer_default = .false.
      if (mprognum_amode(nacc) <= 0) do_aitacc_transfer_default = .false.
   end if

   if ( .not. do_adjust_default ) return

   !  define history fields for number-adjust source-sink for all modes
   do n = 1, ntot_amode
      if (mprognum_amode(n) <= 0) cycle

      do jac = 1, 2
         if (jac == 1) then
            tmpnamea = cnst_name(numptr_amode(n))
         else
            tmpnamea = cnst_name_cw(numptrcw_amode(n))
         end if
         unit = '#/m2/s'
         fieldname = trim(tmpnamea) // '_sfcsiz1'
         long_name = trim(tmpnamea) // ' calcsize number-adjust column source'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnamea) // '_sfcsiz2'
         long_name = trim(tmpnamea) // ' calcsize number-adjust column sink'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname
      end do   ! jac = ...
   end do   ! n = ...

   if ( .not. do_aitacc_transfer_default ) return

   ! check that renaming ipair=1 is aitken-->accum
   ipair = 1
   if ((modefrm_renamexf(ipair) .ne. nait) .or.   &
      (modetoo_renamexf(ipair) .ne. nacc)) then
      write( 6, '(//2a//)' )   &
         '*** modal_aero_calcaersize_init error -- ',   &
         'modefrm/too_renamexf(1) are wrong'
      call endrun( 'modal_aero_calcaersize_init error' )
   end if

   ! define history fields for aitken-accum transfer
   do iq = 1, nspecfrm_renamexf(ipair)

      ! jac=1 does interstitial ("_a"); jac=2 does activated ("_c");
      do jac = 1, 2

         if (jac .eq. 1) then
            lsfrm = lspecfrma_renamexf(iq,ipair)
            lstoo = lspectooa_renamexf(iq,ipair)
         else
            lsfrm = lspecfrmc_renamexf(iq,ipair)
            lstoo = lspectooc_renamexf(iq,ipair)
         end if
         if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

         if (jac .eq. 1) then
            tmpnamea = cnst_name(lsfrm)
            tmpnameb = cnst_name(lstoo)
         else
            tmpnamea = cnst_name_cw(lsfrm)
            tmpnameb = cnst_name_cw(lstoo)
         end if

         unit = 'kg/m2/s'
         if ((tmpnamea(1:3) == 'num') .or. &
            (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
         fieldname = trim(tmpnamea) // '_sfcsiz3'
         long_name = trim(tmpnamea) // ' calcsize aitken-to-accum adjust column tendency'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnameb) // '_sfcsiz3'
         long_name = trim(tmpnameb) // ' calcsize aitken-to-accum adjust column tendency'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnamea) // '_sfcsiz4'
         long_name = trim(tmpnamea) // ' calcsize accum-to-aitken adjust column tendency'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnameb) // '_sfcsiz4'
         long_name = trim(tmpnameb) // ' calcsize accum-to-aitken adjust column tendency'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

      end do   ! jac = ...
   end do   ! iq = ...

#endif

end subroutine modal_aero_calcsize_init

!===============================================================================

subroutine modal_aero_calcsize_sub(state, ptend, deltat, pbuf, aero_props, aero_state, &
   do_adjust_in, do_aitacc_transfer_in)
   ! arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend), target, intent(inout) :: ptend
   real(r8),                    intent(in)    :: deltat
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   class(aerosol_properties), intent(in), target :: aero_props
   class(aerosol_state), intent(in) :: aero_state

   logical, optional :: do_adjust_in
   logical, optional :: do_aitacc_transfer_in

#ifdef MODAL_AERO

   ! local
   logical :: do_adjust
   logical :: do_aitacc_transfer

   integer  :: lchnk, ncol
   integer  :: i, ipair, iq, jac, k, l, lc, lsfrm, lstoo, n
   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: fldcw(:,:)

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname

   ! Work arrays for portable routine interface
   real(r8), allocatable :: q_local(:,:,:)        ! interstitial species (ncol,pver,pcnst)
   real(r8), allocatable :: q_cw_local(:,:,:)     ! cloud-borne species (ncol,pver,pcnst)
   real(r8), allocatable :: dqdt_local(:,:,:)     ! interstitial tendencies (ncol,pver,pcnst)
   real(r8), allocatable :: dqdt_cw_local(:,:,:)  ! cloud-borne tendencies (ncol,pver,pcnst)
   logical,  allocatable :: dotend_local(:)       ! interstitial tendency flags (pcnst)
   logical,  allocatable :: dotend_cw_local(:)    ! cloud-borne tendency flags (pcnst)
   real(r8), allocatable :: qsrflx(:,:,:,:)       ! diagnostic flux (ncol,pcnst,calcsize_nsrflx,2)
   real(r8)              :: qsrflx_pcols(pcols)   ! work array for outfld

   character(len=512)    :: errmsg
   integer               :: errflg

   !-----------------------------------------------------------------------

   if (present(do_adjust_in)) then
      do_adjust = do_adjust_in
   else
      do_adjust = do_adjust_default
   end if

   if (present(do_aitacc_transfer_in)) then
      do_aitacc_transfer = do_aitacc_transfer_in
   else
      do_aitacc_transfer = do_aitacc_transfer_default
   end if

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Get dgnum from pbuf
   call pbuf_get_field(pbuf, dgnum_idx, dgncur_a)

   ! Allocate work arrays
   allocate(q_local(ncol, pver, pcnst))
   allocate(q_cw_local(ncol, pver, pcnst))
   allocate(dqdt_local(ncol, pver, pcnst))
   allocate(dqdt_cw_local(ncol, pver, pcnst))
   allocate(dotend_local(pcnst))
   allocate(dotend_cw_local(pcnst))
   allocate(qsrflx(ncol, pcnst, calcsize_nsrflx, 2))

   ! Copy interstitial species from state%q
   q_local(:,:,:) = state%q(:ncol,:,:)

   ! Gather cloud-borne species from qqcw into flat array
   q_cw_local(:,:,:) = 0.0_r8
   do n = 1, ntot_amode
      do l = 1, nspec_amode(n)
         lc = lmassptrcw_amode(l,n)
         if (lc > 0) then
            fldcw => qqcw_get_field(pbuf, lc, lchnk)
            if (associated(fldcw)) q_cw_local(:ncol,:,lc) = fldcw(:ncol,:)
         end if
      end do
      lc = numptrcw_amode(n)
      if (lc > 0) then
         fldcw => qqcw_get_field(pbuf, lc, lchnk, .true.)
         if (associated(fldcw)) q_cw_local(:ncol,:,lc) = fldcw(:ncol,:)
      end if
   end do

   ! Call portable science routine
   call modal_aero_calcsize_run( &
      ncol               = ncol,                    &
      pver               = pver,                    &
      deltat             = deltat,                  &
      top_lev            = top_lev,                 &
      ntot_amode         = ntot_amode,              &
      nspec_amode        = nspec_amode,             &
      nspec_max          = nspec_max,               &
      dgnum_amode        = dgnum_amode,             &
      dgnumlo_amode      = dgnumlo_amode,           &
      dgnumhi_amode      = dgnumhi_amode,           &
      alnsg_amode        = alnsg_amode,             &
      voltonumb_amode    = voltonumb_amode,         &
      voltonumblo_amode  = voltonumblo_amode,       &
      voltonumbhi_amode  = voltonumbhi_amode,       &
      specdens_amode     = specdens_amode,          &
      mprognum_amode     = mprognum_amode,          &
      modeptr_aitken     = modeptr_aitken,          &
      modeptr_accum      = modeptr_accum,           &
      lmassptr_amode     = lmassptr_amode,          &
      numptr_amode       = numptr_amode,            &
      lmassptrcw_amode   = lmassptrcw_amode,        &
      numptrcw_amode     = numptrcw_amode,          &
      pdel               = state%pdel(:ncol,:),     &
      gravit             = gravit,                  &
      pi                 = pi,                      &
      num_q              = pcnst,                   &
      q                  = q_local,                 &
      q_cw               = q_cw_local,              &
      do_adjust          = do_adjust,               &
      do_aitacc_transfer = do_aitacc_transfer,      &
      npair_renamexf     = npair_renamexf,          &
      nspecfrm_renamexf  = nspecfrm_renamexf,       &
      modefrm_renamexf   = modefrm_renamexf,        &
      modetoo_renamexf   = modetoo_renamexf,        &
      lspecfrma_renamexf = lspecfrma_renamexf,      &
      lspectooa_renamexf = lspectooa_renamexf,      &
      lspecfrmc_renamexf = lspecfrmc_renamexf,      &
      lspectooc_renamexf = lspectooc_renamexf,      &
      dgncur_a           = dgncur_a(:ncol,:,:),     &
      dqdt               = dqdt_local,              &
      dqdt_cw            = dqdt_cw_local,           &
      dotend             = dotend_local,            &
      dotend_cw          = dotend_cw_local,         &
      qsrflx             = qsrflx,                 &
      errmsg             = errmsg,                  &
      errflg             = errflg)

   if (errflg /= 0) then
      call endrun('modal_aero_calcsize_sub: ' // trim(errmsg))
   end if

   ! Apply interstitial tendencies to ptend
   do l = 1, pcnst
      if (dotend_local(l)) then
         ptend%lq(l) = .true.
         ptend%q(:ncol,:,l) = dqdt_local(:,:,l)
      end if
   end do

   ! Apply cloud-borne tendencies to qqcw pbuf fields
   do l = 1, pcnst
      if (dotend_cw_local(l)) then
         fldcw => qqcw_get_field(pbuf, l, lchnk)
         if (associated(fldcw)) then
            do k = top_lev, pver
               do i = 1, ncol
                  fldcw(i,k) = max( 0.0_r8,   &
                     (fldcw(i,k) + dqdt_cw_local(i,k,l)*deltat) )
               end do
            end do
         end if
      end if
   end do

   ! History output
   if (do_adjust) then

      do n = 1, ntot_amode
         if (mprognum_amode(n) <= 0) cycle

         do jac = 1, 2
            if (jac == 1) then
               l = numptr_amode(n)
               tmpnamea = cnst_name(l)
            else
               l = numptrcw_amode(n)
               tmpnamea = cnst_name_cw(l)
            end if

            ! Expand ncol-sized qsrflx to pcols for outfld
            qsrflx_pcols(:) = 0.0_r8
            qsrflx_pcols(:ncol) = qsrflx(:,l,1,jac)
            fieldname = trim(tmpnamea) // '_sfcsiz1'
            call outfld( fieldname, qsrflx_pcols, pcols, lchnk)

            qsrflx_pcols(:) = 0.0_r8
            qsrflx_pcols(:ncol) = qsrflx(:,l,2,jac)
            fieldname = trim(tmpnamea) // '_sfcsiz2'
            call outfld( fieldname, qsrflx_pcols, pcols, lchnk)
         end do   ! jac = ...
      end do   ! n = ...

      if (do_aitacc_transfer) then

         ipair = 1
         do iq = 1, nspecfrm_renamexf(ipair)
            do jac = 1, 2
               if (jac .eq. 1) then
                  lsfrm = lspecfrma_renamexf(iq,ipair)
                  lstoo = lspectooa_renamexf(iq,ipair)
               else
                  lsfrm = lspecfrmc_renamexf(iq,ipair)
                  lstoo = lspectooc_renamexf(iq,ipair)
               end if
               if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

               if (jac .eq. 1) then
                  tmpnamea = cnst_name(lsfrm)
                  tmpnameb = cnst_name(lstoo)
               else
                  tmpnamea = cnst_name_cw(lsfrm)
                  tmpnameb = cnst_name_cw(lstoo)
               end if
               if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

               qsrflx_pcols(:) = 0.0_r8
               qsrflx_pcols(:ncol) = qsrflx(:,lsfrm,3,jac)
               fieldname = trim(tmpnamea) // '_sfcsiz3'
               call outfld( fieldname, qsrflx_pcols, pcols, lchnk)

               qsrflx_pcols(:) = 0.0_r8
               qsrflx_pcols(:ncol) = qsrflx(:,lstoo,3,jac)
               fieldname = trim(tmpnameb) // '_sfcsiz3'
               call outfld( fieldname, qsrflx_pcols, pcols, lchnk)

               qsrflx_pcols(:) = 0.0_r8
               qsrflx_pcols(:ncol) = qsrflx(:,lsfrm,4,jac)
               fieldname = trim(tmpnamea) // '_sfcsiz4'
               call outfld( fieldname, qsrflx_pcols, pcols, lchnk)

               qsrflx_pcols(:) = 0.0_r8
               qsrflx_pcols(:ncol) = qsrflx(:,lstoo,4,jac)
               fieldname = trim(tmpnameb) // '_sfcsiz4'
               call outfld( fieldname, qsrflx_pcols, pcols, lchnk)

            end do   ! jac = ...
         end do   ! iq = ...

      end if   ! do_aitacc_transfer

   end if   ! do_adjust

   call modal_aero_calcdry(state, pbuf, aero_props, aero_state)

   ! Deallocate work arrays
   deallocate(q_local, q_cw_local, dqdt_local, dqdt_cw_local)
   deallocate(dotend_local, dotend_cw_local, qsrflx)

#endif

end subroutine modal_aero_calcsize_sub

subroutine modal_aero_calcsize_diag(state, pbuf, aero_props, aero_state)

   !-----------------------------------------------------------------------
   !
   ! Calculate aerosol size distribution parameters
   !
   ! ***N.B.*** DGNUM for the modes in the climate list are put directly into
   !            the physics buffer.
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), intent(in), target :: state   ! Physics state variables
   type(physics_buffer_desc), pointer :: pbuf(:)      ! physics buffer
   class(aerosol_properties), intent(in), target :: aero_props
   class(aerosol_state), intent(in), target :: aero_state

   ! local
   integer  :: i, k, l1, n
   integer  :: lchnk, ncol
   integer  :: stat
   integer  :: nmodes
   integer  :: nspec

   real(r8), pointer :: dgncur_a(:,:) ! (pcols,pver)


   real(r8), parameter :: third = 1.0_r8/3.0_r8

   real(r8), pointer :: mode_num(:,:) ! mode number mixing ratio
   real(r8), pointer :: specmmr(:,:)  ! specie mmr
   real(r8)          :: specdens      ! specie density

   real(r8) :: dryvol_a(pcols,pver)   ! interstital aerosol dry volume (cm^3/mol_air)

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

   lchnk = state%lchnk
   ncol  = state%ncol

   nmodes = aero_props%nbins()

   do n = 1, nmodes

      call pbuf_get_field(pbuf, dgnum_idx, dgncur_a, start=(/1,1,n/), kount=(/pcols,pver,1/))

      ! get mode properties
      dgnum = aero_props%dgnum(n)
      dgnumhi = aero_props%dgnumhi(n)
      dgnumlo = aero_props%dgnumlo(n)
      sigmag = exp(aero_props%alogsig(n))

      ! get mode number mixing ratio
      call aero_state%get_ambient_num(n, mode_num)

      dgncur_a(:,:) = dgnum
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
                  dgncur_a(i,k) = dgnxx
               else if (num_a >= drv_a*v2nyy) then
                  dgncur_a(i,k) = dgnyy
               else
                  dgncur_a(i,k) = (drv_a/(dumfac*num_a))**third
               end if
            end if

         end do
      end do

   end do ! nmodes

   call modal_aero_calcdry(state, pbuf, aero_props, aero_state)

end subroutine modal_aero_calcsize_diag

subroutine modal_aero_calcdry(state, pbuf, aero_props, aero_state)

   type(physics_state), target, intent(in)    :: state       ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer
   class(aerosol_properties), intent(in), target :: aero_props
   class(aerosol_state), intent(in), target :: aero_state

   real(r8), pointer :: hygro(:,:,:)
   real(r8), pointer :: dryvol(:,:,:)
   real(r8), pointer :: dryrad(:,:,:)
   real(r8), pointer :: drymass(:,:,:)
   real(r8), pointer :: so4dryvol(:,:,:)
   real(r8), pointer :: naer(:,:,:)
   real(r8), pointer :: dgncur_a(:,:,:)

   integer :: ncol

   character(len=512) :: errmsg_local
   integer            :: errflg_local

   ncol = state%ncol

   call pbuf_get_field(pbuf, dgnum_idx,     dgncur_a)
   call pbuf_get_field(pbuf, hygro_idx,     hygro)
   call pbuf_get_field(pbuf, dryvol_idx,    dryvol)
   call pbuf_get_field(pbuf, dryrad_idx,    dryrad)
   call pbuf_get_field(pbuf, drymass_idx,   drymass)
   call pbuf_get_field(pbuf, so4dryvol_idx, so4dryvol)
   call pbuf_get_field(pbuf, naer_idx,      naer)

   hygro(:,:,:)     = 0._r8
   dryvol(:,:,:)    = 0._r8
   dryrad(:,:,:)    = 0._r8
   drymass(:,:,:)   = 0._r8
   so4dryvol(:,:,:) = 0._r8
   naer(:,:,:)      = 0._r8

   ! call portable subroutine:
   call modal_aero_calcdry_run( &
      aero_props       = aero_props,            &
      aero_state       = aero_state,            &
      ncol             = ncol,                  &
      pver             = pver,                  &
      top_lev          = top_lev,               &
      do_strat_sulfate = modal_strat_sulfate,   &
      pi               = pi,                    &
      dgncur_a         = dgncur_a(:ncol,:,:),   &
      hygro            = hygro(:ncol,:,:),       &
      dryvol           = dryvol(:ncol,:,:),      &
      dryrad           = dryrad(:ncol,:,:),      &
      drymass          = drymass(:ncol,:,:),     &
      so4dryvol        = so4dryvol(:ncol,:,:),   &
      naer             = naer(:ncol,:,:),        &
      errmsg           = errmsg_local,          &
      errflg           = errflg_local)
   if (errflg_local /= 0) then
      call endrun('modal_aero_calcdry: ' // trim(errmsg_local))
   end if

end subroutine modal_aero_calcdry

end module modal_aero_calcsize_cam
