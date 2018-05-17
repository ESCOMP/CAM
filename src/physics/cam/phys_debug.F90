module phys_debug

! This module contains subroutines that are intended to print diagnostic
! information about a single column to the log file.  The module uses the
! phys_debug_util module to locate the column using local indices of the
! physics decomposition.  The reason to encapsulate the debug print
! statements into a subroutine is to minimize the impact of the debugging
! code on the development code.  This makes it easier to maintain debug
! code on a branch which can be updated to the head of the trunk when
! convenient.  Minimizing the footprint of debug code helps to minimize the
! chances of conflicts that can occur while updating between the debug code
! and changes that have occurred on the trunk.

! This module also contains some generic code to write a collection of variables
! from the physics state object to the history file.


use shr_kind_mod,    only: r8 => shr_kind_r8
use phys_debug_util, only: phys_debug_col
use ppgrid,          only: pcols, pver, pverp
use physics_types,   only: physics_state, physics_ptend
use camsrfexch,      only: cam_in_t     
use cam_logfile,     only: iulog
use cam_history,     only: addfld, add_default, outfld
use constituents,    only: pcnst, cnst_get_ind, cnst_name
use rad_constituents,only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx

implicit none
private
save

! Generic routines for writing physics state to history
public :: &
   phys_debug_state_init, &
   phys_debug_state_out


! The following routines are intended as examples for printing various diagnostic
! information from a single column
public :: &
   phys_debug_vdiff1,   &
   phys_debug_shallow1, &
   phys_debug_strat1,   &
   phys_debug_srf1,     &
   phys_debug_hbdiff1,  &
   phys_debug_flux1,    &
   phys_debug_flux2

! constituent indices
integer :: ixcldliq, ixcldice, ixnumliq, ixnumice
! logical array to identify constituents that are mode number concentrations
logical :: cnst_is_mam_num(pcnst)
! logical array to identify constituents that are mode specie mass mixing ratios
logical :: cnst_is_mam_mmr(pcnst)


!===================================================================================================
contains
!===================================================================================================

subroutine phys_debug_state_init(tags)

   ! Write selected state variables from locations corresponding to the 'tags'.
   ! This method does the addfld/add_default calls.

   ! ** Note ** This init needs to be called after the rad_constituent module has finished
   !            its init.


   character(len=8), dimension(:), intent(in) :: tags

   ! local variables
   integer :: i, icnst, j, nmodes, nspec, ntags
   character(len=5) :: unit
   !-----------------------------------------------------------------------------------------------

   ! Initialize some constituent indices and names
   call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
   call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
   call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)
   call cnst_get_ind('NUMICE', ixnumice, abort=.false.)

   ! Set arrays to identify the modal aerosol constituents
   cnst_is_mam_num = .false.
   cnst_is_mam_mmr = .false.
   call rad_cnst_get_info(0, nmodes=nmodes)
   do i = 1, nmodes
      call rad_cnst_get_mode_num_idx(i, icnst)
      cnst_is_mam_num(icnst) = .true.
      call rad_cnst_get_info(0, i, nspec=nspec)
      do j = 1, nspec
         call rad_cnst_get_mam_mmr_idx(i, j, icnst)
         cnst_is_mam_mmr(icnst) = .true.
      end do
   end do

   ! For each supplied tag, issue addfld/add_default calls for the state variables
   ! with the tag name appended.
   ntags = size(tags)
   do i = 1, ntags

      call addfld('T_'//trim(tags(i)), (/ 'lev' /), 'A', 'K', 'T at tag '//trim(tags(i)))
      call add_default('T_'//trim(tags(i)), 1, ' ')

      call addfld('U_'//trim(tags(i)), (/ 'lev' /), 'A', 'm/s', 'U at tag '//trim(tags(i)))
      call add_default('U_'//trim(tags(i)), 1, ' ')

      call addfld('V_'//trim(tags(i)), (/ 'lev' /), 'A', 'm/s', 'V at tag '//trim(tags(i)))
      call add_default('V_'//trim(tags(i)), 1, ' ')

      call addfld('QV_'//trim(tags(i)), (/ 'lev' /), 'A', 'kg/kg', 'QV at tag '//trim(tags(i)))
      call add_default('QV_'//trim(tags(i)), 1, ' ')

      if (ixcldliq > 0) then
         call addfld('QL_'//trim(tags(i)), (/ 'lev' /), 'A', 'kg/kg', 'QL at tag '//trim(tags(i)))
         call add_default('QL_'//trim(tags(i)), 1, ' ')
      end if

      if (ixcldice > 0) then
         call addfld('QI_'//trim(tags(i)), (/ 'lev' /), 'A', 'kg/kg', 'QI at tag '//trim(tags(i)))
         call add_default('QI_'//trim(tags(i)), 1, ' ')
      end if

      if (ixnumliq > 0) then
         call addfld('NL_'//trim(tags(i)), (/ 'lev' /), 'A', '1/kg', 'NL at tag '//trim(tags(i)))
         call add_default('NL_'//trim(tags(i)), 1, ' ')
      end if

      if (ixnumice > 0) then
         call addfld('NI_'//trim(tags(i)), (/ 'lev' /), 'A', '1/kg', 'NI at tag '//trim(tags(i)))
         call add_default('NI_'//trim(tags(i)), 1, ' ')
      end if

      do j = 1, pcnst

         if (cnst_is_mam_num(j) .or. cnst_is_mam_mmr(j)) then

            unit = '1/kg'
            if (cnst_is_mam_mmr(j)) unit = 'kg/kg'

            call addfld(trim(cnst_name(j))//'_'//trim(tags(i)), (/ 'lev' /), &
               'A', trim(unit), trim(cnst_name(j))//' at tag '//trim(tags(i)))
            call add_default(trim(cnst_name(j))//'_'//trim(tags(i)), 1, ' ')

         end if
      end do

   end do

end subroutine phys_debug_state_init

!===================================================================================================

subroutine phys_debug_state_out(tag, state)

   ! Dump selected state variables at location corresponding to given tag.

   character(len=*),    intent(in) :: tag
   type(physics_state), intent(in) :: state 

   ! Local variables
   integer :: lchnk, j
   !-----------------------------------------------------------------------------------------------

   lchnk = state%lchnk

   call outfld('T_'//trim(tag),  state%t,        pcols, lchnk)
   call outfld('U_'//trim(tag),  state%u,        pcols, lchnk)
   call outfld('V_'//trim(tag),  state%v,        pcols, lchnk)
   call outfld('QV_'//trim(tag), state%q(:,:,1), pcols, lchnk)
   if (ixcldliq > 0) call outfld('QL_'//trim(tag), state%q(:,:,ixcldliq), pcols, lchnk)
   if (ixcldice > 0) call outfld('QI_'//trim(tag), state%q(:,:,ixcldice), pcols, lchnk)
   if (ixnumliq > 0) call outfld('NL_'//trim(tag), state%q(:,:,ixnumliq), pcols, lchnk)
   if (ixnumice > 0) call outfld('NI_'//trim(tag), state%q(:,:,ixnumice), pcols, lchnk)

   do j = 1, pcnst
      if (cnst_is_mam_num(j) .or. cnst_is_mam_mmr(j)) then
         call outfld(trim(cnst_name(j))//'_'//trim(tag), state%q(:,:,j), pcols, lchnk)
      end if
   end do

end subroutine phys_debug_state_out


!===================================================================================================

subroutine phys_debug_vdiff1(state, kvm, tautotx, tautoty)

   type(physics_state), intent(in) :: state             ! Physics state variables
   real(r8),            intent(in) :: kvm(pcols,pverp)  ! eddy diffusivity for momentum [m2/s]
   real(r8),            intent(in) :: tautotx(pcols)    ! u component of total surface stress
   real(r8),            intent(in) :: tautoty(pcols)    ! v component of total surface stress

   integer :: icol
   !-----------------------------------------------------------------------------

   icol = phys_debug_col(state%lchnk) 
   if (icol > 0) then
      write(iulog,*) ' vert_diff: kvm, u ', kvm(icol,pver), kvm(icol,pver-1), &
         state%u(icol, pver), state%u(icol,pver-1), &
         state%v(icol, pver), state%v(icol,pver-1)
      write(iulog,*) ' vert_diff: taux, tauy ', tautotx(icol), tautoty(icol)
   endif

end subroutine phys_debug_vdiff1

!================================================================================

subroutine phys_debug_shallow1(state, ptend, nstep, prec_cmf, rliq2, ztodt, kmx)

   use constituents, only: cnst_get_ind

   type(physics_state), intent(in) :: state             ! Physics state variables
   type(physics_ptend), intent(in) :: ptend             ! Physics process tendencies
   integer,             intent(in) :: nstep
   real(r8),            intent(in) :: prec_cmf(pcols)   ! total precipitation from Hack convection
   real(r8),            intent(in) :: rliq2(pcols)      ! vertical integral of liquid from shallow scheme
   real(r8),            intent(in) :: ztodt             ! physics time step
   integer,            intent(out) :: kmx

   integer  :: icol, k, kmn, ixcldliq
   real(r8) :: qtmx, qtmn
   !-----------------------------------------------------------------------------

   icol = phys_debug_col(state%lchnk) 
   if (icol > 0) then
      qtmx = 0._r8
      kmx = 0
      qtmn = 0._r8
      kmn = 0
      do k = 1,pver
         !        write (iulog,*) 'aaa ', k, ptend%q(icol,k,1), ptend%q(icol,k,2)
         if (ptend%q(icol,k,1).ge.qtmx) then
            kmx = k
            qtmx = ptend%q(icol,k,1)
         endif
         if (ptend%q(icol,k,1).le.qtmn) then
            kmn = k
            qtmn = ptend%q(icol,k,1)
         endif
      end do
      k = kmx
66    format ('tphysbc, aft shallow:', 4i5, 6f9.4) 
      call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
      write (iulog,66)  nstep, icol, &
           kmx, kmn, &
           prec_cmf(icol)*8.64e7_r8, rliq2(icol)*8.64e7_r8,  &
           qtmx*8.64e7_r8, qtmn*8.64e7_r8, &
           (state%q(icol,k,1)+ptend%q(icol,k,1)*ztodt)*1.e3_r8, &
           (state%q(icol,k,ixcldliq)+ptend%q(icol,k,ixcldliq)*ztodt)*1.e3_r8

   endif

end subroutine phys_debug_shallow1

!================================================================================

subroutine phys_debug_strat1(state, ptend, nstep, prec_str, rliq, ztodt, kmx)

   use constituents, only: cnst_get_ind

   type(physics_state), intent(in) :: state             ! Physics state variables
   type(physics_ptend), intent(in) :: ptend             ! Physics process tendencies
   integer,             intent(in) :: nstep
   real(r8),            intent(in) :: prec_str(pcols)   ! sfc flux of precip from stratiform (m/s)
   real(r8),            intent(in) :: rliq(pcols)       ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8),            intent(in) :: ztodt             ! physics time step
   integer,             intent(in) :: kmx

   integer  :: icol, k, ixcldliq
   !-----------------------------------------------------------------------------

   icol = phys_debug_col(state%lchnk) 
   if (icol > 0) then

      k = kmx
67    format ('tphysbc, aft strat:', i5, 6f9.4) 
      call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
      write (iulog,67)  nstep, prec_str(icol)*8.64e7_r8, rliq(icol)*8.64e7_r8,  &
           (ptend%q(icol,k,1)*ztodt)*1.e3_r8, &
           (ptend%q(icol,k,ixcldliq)*ztodt)*1.e3_r8, &
           (state%q(icol,k,1)+ptend%q(icol,k,1)*ztodt)*1.e3_r8, &
           (state%q(icol,k,ixcldliq)+ptend%q(icol,k,ixcldliq)*ztodt)*1.e3_r8

   endif

end subroutine phys_debug_strat1

!================================================================================

subroutine phys_debug_srf1(lchnk, cam_in)

   integer,             intent(in) :: lchnk             ! local chunk index
   type(cam_in_t),      intent(in) :: cam_in            ! CAM's import state

   integer :: icol
   !-----------------------------------------------------------------------------

   icol = phys_debug_col(lchnk) 
   if (icol > 0) then

      write(iulog,*) 'bot tphysbc: cam_in%tref', cam_in%tref(icol)

   endif

end subroutine phys_debug_srf1

!================================================================================

subroutine phys_debug_hbdiff1(lchnk, pblh, zl, zh)

   integer,  intent(in) :: lchnk             ! local chunk index
   real(r8), intent(in) :: pblh(pcols)       ! boundary-layer height [m]
   real(r8), intent(in) :: zl(pcols)         ! zmzp / Obukhov length
   real(r8), intent(in) :: zh(pcols)         ! zmzp / pblh

   integer :: icol
   !-----------------------------------------------------------------------------

   icol = phys_debug_col(lchnk) 
   if (icol > 0) then

      write(iulog,*) ' austach_pbl, pblh, zl, zh: ', pblh(icol), zl(icol), zh(icol)

   endif

end subroutine phys_debug_hbdiff1

!================================================================================

subroutine phys_debug_flux1(lchnk, srfflx, lhflx, shflx, taux, tauy, qflx, &
                            lhflx_res, shflx_res, taux_res, tauy_res, qflx_res)

   integer,            intent(in) :: lchnk      ! local chunk index
   type(cam_in_t),     intent(in) :: srfflx     ! cam import state
   real(r8),           intent(in) :: lhflx(:)   ! latent heat flux
   real(r8),           intent(in) :: shflx(:)   ! sensible heat flux
   real(r8),           intent(in) :: taux(:)    ! x momentum flux
   real(r8),           intent(in) :: tauy(:)    ! y momentum flux
   real(r8),           intent(in) :: qflx(:)    ! water vapor heat flux
   real(r8),           intent(in) :: lhflx_res(:)   ! latent heat flux
   real(r8),           intent(in) :: shflx_res(:)   ! sensible heat flux
   real(r8),           intent(in) :: taux_res(:)    ! x momentum flux
   real(r8),           intent(in) :: tauy_res(:)    ! y momentum flux
   real(r8),           intent(in) :: qflx_res(:)    ! water vapor heat flux

   integer :: icol
   !-----------------------------------------------------------------------------

   icol = phys_debug_col(lchnk) 
   if (icol > 0) then

      write(iulog,*) ' b flux_tweak called, lhflx, oldlhflx ', &
          srfflx%lhf(icol), lhflx(icol)
      write(iulog,*) ' sfmodel fluxes lhf, shf, taux, tauy, q: ', &
          srfflx%lhf(icol), &
          srfflx%shf(icol), &
          srfflx%wsx(icol), &
          srfflx%wsy(icol), &
          srfflx%cflx(icol,1)
      write(iulog,*) ' last fluxes used lhf, shf, taux, tauy, q: ', &
          lhflx(icol), &
          shflx(icol), &
          taux(icol), &
          tauy(icol), &
          qflx(icol)
      write(iulog,*) ' current residuals lhf, shf, taux, tauy, q: ', &
          lhflx_res(icol), &
          shflx_res(icol), &
          taux_res(icol), &
          tauy_res(icol), &
          qflx_res(icol)

   endif

end subroutine phys_debug_flux1

!================================================================================

subroutine phys_debug_flux2(lchnk, srfflx, lhflx, &
                            lhflx_res, shflx_res, taux_res, tauy_res, qflx_res)

   integer,            intent(in) :: lchnk      ! local chunk index
   type(cam_in_t),     intent(in) :: srfflx     ! cam import state
   real(r8),           intent(in) :: lhflx(:)   ! latent heat flux
   real(r8),           intent(in) :: lhflx_res(:)   ! latent heat flux
   real(r8),           intent(in) :: shflx_res(:)   ! sensible heat flux
   real(r8),           intent(in) :: taux_res(:)    ! x momentum flux
   real(r8),           intent(in) :: tauy_res(:)    ! y momentum flux
   real(r8),           intent(in) :: qflx_res(:)    ! water vapor heat flux

   integer :: icol
   !-----------------------------------------------------------------------------

   icol = phys_debug_col(lchnk) 
   if (icol > 0) then

      write(iulog,*) ' a flux_tweak called, lhflx, oldlhflx ', &
          srfflx%lhf(icol), lhflx(icol)
      write (iulog,66) ' residual fractions lhf, shf, taux, tauy, q ', &
          lhflx_res(icol)/srfflx%lhf(icol), &
          shflx_res(icol)/srfflx%shf(icol), &
          taux_res(icol)/srfflx%wsx(icol), &
          tauy_res(icol)/srfflx%wsy(icol), &
          qflx_res(icol)/srfflx%cflx(icol,1)
      write(iulog,66) ' residual lhf, shf, taux, tauy, q: ', &
          lhflx_res(icol), &
          shflx_res(icol), &
          taux_res(icol), &
          tauy_res(icol), &
          qflx_res(icol)
      write(iulog,*) ' used fluxes lhf, shf, taux, tauy, q: ', &
          srfflx%lhf(icol), &
          srfflx%shf(icol), &
          srfflx%wsx(icol), &
          srfflx%wsy(icol), &
          srfflx%cflx(icol,1)
66 format (a, 1p, 5e15.5)

   endif

end subroutine phys_debug_flux2

end module phys_debug
