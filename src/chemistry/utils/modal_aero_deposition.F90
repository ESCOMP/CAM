module modal_aero_deposition

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from modal components of wet and dry 
! deposition at the surface into the fields passed to the coupler.
!
! *** N.B. *** Currently only a simple scheme for the 3-mode version
!              of MAM has been implemented.
!
! Revision history:
! Feb 2009  M. Flanner, B. Eaton   Original version for trop_mam3.
! Jul 2011  F Vitt -- made avaliable to be used in a prescribed modal aerosol mode (no prognostic MAM)
! Mar 2012  F Vitt -- made changes for to prevent abort when 7-mode aeroslol model is used
!                     some of the needed consituents do not exist in 7-mode so bin_fluxes will be false
! May 2014  F Vitt -- included contributions from MAM4 aerosols and added soa_a2 to the ocphiwet fluxes
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use camsrfexch,       only: cam_out_t     
use constituents,     only: cnst_get_ind, pcnst
use cam_abortutils,   only: endrun
use rad_constituents, only: rad_cnst_get_info

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep

! Private module data

logical :: initialized = .false.
integer :: bcphi_ndx( pcnst ) = -1
integer :: bcpho_ndx( pcnst ) = -1
integer :: ocphi_ndx( pcnst ) = -1
integer :: ocpho_ndx( pcnst ) = -1
integer :: crse_dust_ndx( pcnst ) = -1
integer :: fine_dust_ndx( pcnst ) = -1
integer :: bcphi_cnt = 0
integer :: ocphi_cnt = 0
integer :: bcpho_cnt = 0
integer :: ocpho_cnt = 0
integer :: crse_dust_cnt = 0
integer :: fine_dust_cnt = 0

!==============================================================================
contains
!==============================================================================

subroutine modal_aero_deposition_init( bcphi_indices, bcpho_indices, ocphi_indices, &
                                ocpho_indices, fine_dust_indices, crse_dust_indices )

  ! set aerosol indices for re-mapping surface deposition fluxes:
  ! *_a1 = accumulation mode
  ! *_a2 = aitken mode
  ! *_a3 = coarse mode
  
  ! can be initialized with user specified indices
  ! if called from aerodep_flx module (for prescribed modal aerosol fluxes) then these indices are specified
  integer, optional, intent(in) :: bcphi_indices(:)     ! hydrophilic black carbon
  integer, optional, intent(in) :: bcpho_indices(:)     ! hydrophobic black carbon
  integer, optional, intent(in) :: ocphi_indices(:)     ! hydrophilic organic carbon
  integer, optional, intent(in) :: ocpho_indices(:)     ! hydrophobic organic carbon 
  integer, optional, intent(in) :: fine_dust_indices(:) ! fine dust
  integer, optional, intent(in) :: crse_dust_indices(:) ! coarse dust

  ! local vars
  integer :: i, pcnt, scnt

  character(len=16), parameter :: fine_dust_modes(2) =  (/ 'accum           ', 'fine_dust       '/)
  character(len=16), parameter :: crse_dust_modes(2) =  (/ 'coarse          ', 'coarse_dust     '/)
  character(len=16), parameter :: hydrophilic_carbon_modes(1) = (/'accum           '/)
  character(len=16), parameter :: hydrophobic_carbon_modes(3) = (/'aitken          ',  'coarse          ', 'primary_carbon  '/)

  ! if already initialized abort the run
  if (initialized) then
     call endrun('modal_aero_deposition is already initialized')
  endif

  if (present(bcphi_indices)) then
     bcphi_cnt = size(bcphi_indices)
     bcphi_ndx(1:bcphi_cnt) = bcphi_indices (1:bcphi_cnt)
  else
     call get_indices( type='black-c', modes=hydrophilic_carbon_modes, indices=bcphi_ndx, count=bcphi_cnt )
  endif
  if (present(bcpho_indices)) then
     bcpho_cnt = size(bcpho_indices)
     bcpho_ndx(1:bcpho_cnt) = bcpho_indices (1:bcpho_cnt)
  else
     call get_indices( type='black-c', modes=hydrophobic_carbon_modes, indices=bcpho_ndx, count=bcpho_cnt )
  endif

  if (present(ocphi_indices)) then
     ocphi_cnt = size(ocphi_indices)
     ocphi_ndx(1:ocphi_cnt) = ocphi_indices (1:ocphi_cnt)
  else
     call get_indices( type='s-organic', modes=hydrophilic_carbon_modes, indices=ocphi_ndx, count=pcnt )
     call get_indices( type='p-organic', modes=hydrophilic_carbon_modes, indices=ocphi_ndx(pcnt+1:), count=scnt )
     ocphi_cnt = pcnt+scnt
  endif
  if (present(ocpho_indices)) then
     ocpho_cnt = size(ocpho_indices)
     ocpho_ndx(1:ocpho_cnt) = ocpho_indices (1:ocpho_cnt)
  else
     call get_indices( type='s-organic', modes=hydrophobic_carbon_modes, indices=ocpho_ndx, count=pcnt )
     call get_indices( type='p-organic', modes=hydrophobic_carbon_modes, indices=ocpho_ndx(pcnt+1:), count=scnt )
     ocpho_cnt = pcnt+scnt
  endif

  if (present(fine_dust_indices)) then
     fine_dust_cnt = size(fine_dust_indices)
     fine_dust_ndx(1:fine_dust_cnt) = fine_dust_indices(1:fine_dust_cnt)
  else
     call get_indices( type='dust', modes=fine_dust_modes, indices=fine_dust_ndx, count=fine_dust_cnt )
  endif
  if (present(crse_dust_indices)) then
     crse_dust_cnt = size(crse_dust_indices)
     crse_dust_ndx(1:crse_dust_cnt) = crse_dust_indices(1:crse_dust_cnt)
  else
     call get_indices( type='dust', modes=crse_dust_modes, indices=crse_dust_ndx, count=crse_dust_cnt )
  endif

  initialized = .true.

end subroutine modal_aero_deposition_init

!==============================================================================
subroutine set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec, idx
   integer :: ncol                      ! number of columns

   real(r8) :: bcphiwet_sum, ocphiwet_sum
   !----------------------------------------------------------------------------

  if (.not.initialized) call endrun('set_srf_wetdep: modal_aero_deposition has not been initialized')

   ncol = cam_out%ncol

   cam_out%bcphiwet(:) = 0._r8
   cam_out%ocphiwet(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      do ispec=1,bcphi_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcphi_ndx(ispec))+aerdepwetcw(i,bcphi_ndx(ispec)))
      enddo
      do ispec=1,bcpho_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcpho_ndx(ispec))+aerdepwetcw(i,bcpho_ndx(ispec)))
      enddo

      ! organic carbon fluxes
      do ispec=1,ocphi_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocphi_ndx(ispec))+aerdepwetcw(i,ocphi_ndx(ispec)))
      enddo
      do ispec=1,ocpho_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocpho_ndx(ispec))+aerdepwetcw(i,ocpho_ndx(ispec)))
      enddo

      ! dust fluxes
      cam_out%dstwet1(i) = 0._r8
      cam_out%dstwet2(i) = 0._r8
      cam_out%dstwet3(i) = 0._r8
      cam_out%dstwet4(i) = 0._r8

      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      do ispec=1,fine_dust_cnt
         cam_out%dstwet1(i) = cam_out%dstwet1(i) &
                            -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
      enddo

      !  Assign all coarse-mode dust to bulk size bin 3:
      do ispec=1,crse_dust_cnt
         cam_out%dstwet3(i) = cam_out%dstwet3(i) &
                            -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
      enddo

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8
   enddo

end subroutine set_srf_wetdep

!==============================================================================

subroutine set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)

! Set surface dry deposition fluxes passed to coupler.
   
   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec, idx
   integer :: ncol                      ! number of columns
   real(r8):: bcphidry_sum, ocphidry_sum, ocphodry_sum
   !----------------------------------------------------------------------------

   if (.not.initialized) call endrun('set_srf_drydep: modal_aero_deposition has not been initialized')

   ncol = cam_out%ncol

   cam_out%bcphidry(:) = 0._r8
   cam_out%bcphodry(:) = 0._r8
   cam_out%ocphidry(:) = 0._r8
   cam_out%ocphodry(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      do ispec=1,bcphi_cnt
         cam_out%bcphidry(i) = cam_out%bcphidry(i) &
                             + (aerdepdryis(i,bcphi_ndx(ispec))+aerdepdrycw(i,bcphi_ndx(ispec)))
      enddo
      do ispec=1,bcpho_cnt
         cam_out%bcphodry(i) = cam_out%bcphodry(i) &
                             + (aerdepdryis(i,bcpho_ndx(ispec))+aerdepdrycw(i,bcpho_ndx(ispec)))
      enddo

      ! organic carbon fluxes
      do ispec=1,ocphi_cnt
         cam_out%ocphidry(i) = cam_out%ocphidry(i) &
                             + (aerdepdryis(i,ocphi_ndx(ispec))+aerdepdrycw(i,ocphi_ndx(ispec)))
      enddo
      do ispec=1,ocpho_cnt
         cam_out%ocphodry(i) = cam_out%ocphodry(i) &
                             + (aerdepdryis(i,ocpho_ndx(ispec))+aerdepdrycw(i,ocpho_ndx(ispec)))
      enddo

      ! dust fluxes
      cam_out%dstdry1(i) = 0._r8
      cam_out%dstdry2(i) = 0._r8
      cam_out%dstdry3(i) = 0._r8
      cam_out%dstdry4(i) = 0._r8
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      do ispec=1,fine_dust_cnt
         cam_out%dstdry1(i) = cam_out%dstdry1(i) &
                            + (aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
      enddo
      !  Assign all coarse-mode dust to bulk size bin 3:
      do ispec=1,crse_dust_cnt
         cam_out%dstdry3(i) = cam_out%dstdry3(i) &
                            + (aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
      enddo

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(i) .lt. 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%bcphodry(i) .lt. 0._r8) cam_out%bcphodry(i) = 0._r8
      if (cam_out%ocphidry(i) .lt. 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%ocphodry(i) .lt. 0._r8) cam_out%ocphodry(i) = 0._r8
      if (cam_out%dstdry1(i)  .lt. 0._r8) cam_out%dstdry1(i)  = 0._r8
      if (cam_out%dstdry3(i)  .lt. 0._r8) cam_out%dstdry3(i)  = 0._r8
   enddo

end subroutine set_srf_drydep

!==============================================================================
subroutine get_indices( type, modes, indices, count )

  character(len=*), intent(in) :: type
  character(len=*), intent(in) :: modes(:)
  integer, intent(out) :: indices(:)
  integer, intent(out) :: count

  integer :: l, n, ndx, nmodes, nspec
  character(len=32) :: spec_type, spec_name, mode_type

  call rad_cnst_get_info(0, nmodes=nmodes)

  count = 0
  indices(:) = -1

  if (nmodes==7) return ! historically turned off for mam7

  do n = 1, nmodes

     call rad_cnst_get_info(0, n, mode_type=mode_type, nspec=nspec)

     if ( any(modes==trim(mode_type)) ) then

        do l = 1,nspec
           call rad_cnst_get_info(0, n, l, spec_type=spec_type, spec_name=spec_name)
           call cnst_get_ind(spec_name, ndx, abort=.false.)
           if (ndx>0) then
              if (trim(spec_type) == trim(type)) then
                 count = count+1
                 indices(count) = ndx
              endif
           endif
        enddo

     endif

  enddo

end subroutine get_indices
!==============================================================================

end module modal_aero_deposition
