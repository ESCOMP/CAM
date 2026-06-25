module modal_aero_wateruptake_cam

!  CAM wrapper for modal_aero_wateruptake.
!  Handles pbuf registration, initialization, history output,
!  state/pbuf/aero_props marshaling, and calls the portable
!  science routines.

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o, rair
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use aerosol_properties_mod, only: aerosol_properties
use aerosol_state_mod, only: aerosol_state
use cam_history,      only: addfld, add_default, outfld, horiz_only
use cam_logfile,      only: iulog
use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use phys_control,     only: phys_getopts
use cam_abortutils,   only: endrun

use modal_aero_wateruptake, only: modal_aero_wateruptake_init, &
                                  modal_aero_wateruptake_sub
use modal_aero_data,          only: modal_strat_sulfate

implicit none
private
save

public :: &
   modal_aero_wateruptake_cam_init, &
   modal_aero_wateruptake_dr,   &
   modal_aero_wateruptake_reg

! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: sulfeq_idx     = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0
integer :: hygro_idx      = 0
integer :: dryvol_idx     = 0
integer :: dryrad_idx     = 0
integer :: drymass_idx    = 0
integer :: so4dryvol_idx  = 0
integer :: naer_idx       = 0


!===============================================================================
contains
!===============================================================================

subroutine modal_aero_wateruptake_reg()

  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use radiative_aerosol, only: rad_aer_get_info

   integer :: nmodes

   call rad_aer_get_info(0, nmodes=nmodes)
   call pbuf_add_field('DGNUMWET',   'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_idx)
   call pbuf_add_field('WETDENS_AP', 'physpkg', dtype_r8, (/pcols, pver, nmodes/), wetdens_ap_idx)

   ! 1st order rate for direct conversion of strat. cloud water to precip (1/s)
   call pbuf_add_field('QAERWAT',    'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_idx)

   if (modal_strat_sulfate) then
      call pbuf_add_field('MAMH2SO4EQ', 'global',  dtype_r8, (/pcols, pver, nmodes/), sulfeq_idx)
   end if


end subroutine modal_aero_wateruptake_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_wateruptake_cam_init(pbuf2d)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field
   use infnan,       only : nan, assignment(=)
   use radiative_aerosol, only: rad_aer_get_info
   use modal_aero_wateruptake, only: modal_aero_wateruptake_diag
   use modal_aerosol_state_mod, only: modal_aerosol_state_register_water_uptake_diag

   use shr_const_mod, only: shr_const_pi

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   real(r8) :: real_nan

   integer :: m, nmodes
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies

   character(len=3) :: trnum       ! used to hold mode number (as characters)

   character(len=512) :: errmsg_local
   integer            :: errflg_local
   !----------------------------------------------------------------------------

   real_nan = nan

   cld_idx        = pbuf_get_index('CLD')
   dgnum_idx      = pbuf_get_index('DGNUM')

   hygro_idx      = pbuf_get_index('HYGRO')
   dryvol_idx     = pbuf_get_index('DRYVOL')
   dryrad_idx     = pbuf_get_index('DRYRAD')
   drymass_idx    = pbuf_get_index('DRYMASS')
   so4dryvol_idx  = pbuf_get_index('SO4DRYVOL')
   naer_idx       = pbuf_get_index('NAER')

   ! assume for now that will compute wateruptake for climate list modes only

   call rad_aer_get_info(0, nmodes=nmodes)

   do m = 1, nmodes
      write(trnum, '(i3.3)') m
      call addfld('dgnd_a'//trnum(2:3), (/ 'lev' /), 'A', 'm', &
         'dry dgnum, interstitial, mode '//trnum(2:3))
      call addfld('dgnw_a'//trnum(2:3), (/ 'lev' /), 'A', 'm', &
         'wet dgnum, interstitial, mode '//trnum(2:3))
      call addfld('wat_a'//trnum(3:3),  (/ 'lev' /), 'A', 'm', &
         'aerosol water, interstitial, mode '//trnum(2:3))

      ! determine default variables
      call phys_getopts(history_aerosol_out = history_aerosol)

      if (history_aerosol) then
         call add_default('dgnd_a'//trnum(2:3), 1, ' ')
         call add_default('dgnw_a'//trnum(2:3), 1, ' ')
         call add_default('wat_a'//trnum(3:3),  1, ' ')
      endif

   end do

   call addfld('PM25',     (/ 'lev' /), 'A', 'kg/m3', 'PM2.5 mass concentration')
   call addfld('PM25_SRF', horiz_only,  'A', 'kg/m3', 'surface PM2.5 mass concentration')
   ! dmleung added a few more below, 20 Nov 2023
   call addfld('PM25_MMR',     (/ 'lev' /), 'A', 'kg/kg', 'PM2.5 mass mixing ratio')
   call addfld('PM1_SRF',      horiz_only,  'A', 'kg/m3', 'surface PM1 mass concentration')
   call addfld('PM1_MMR',      (/ 'lev' /), 'A', 'kg/kg', 'PM1 mass mixing ratio')
   call addfld('PM10_SRF',      horiz_only,  'A', 'kg/m3', 'surface PM10 mass concentration')
   call addfld('PM10_MMR',     (/ 'lev' /), 'A', 'kg/kg', 'PM10 mass mixing ratio')
   call addfld('PMTOT_MMR',     (/ 'lev' /), 'A', 'kg/kg', 'total PM mass mixing ratio')
   call addfld('RHO_AIR',      (/ 'lev' /), 'A', 'kg/m3', 'air density')  ! I know RHO_CLUBB exists. Does this exist?

   call add_default('RHO_AIR', 1, ' ')
   call add_default('PM25_SRF', 1, ' ')
   call add_default('PM25_MMR',  1, ' ')
   call add_default('PM10_MMR',  1, ' ')
   ! dmleung --

   if (is_first_step()) then
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnumwet_idx, 0.0_r8)
      if (modal_strat_sulfate) then
      ! initialize fields in physics buffer to NaN (not a number)
      ! so model will crash if used before initialization
         call pbuf_set_field(pbuf2d, sulfeq_idx, real_nan)
      endif
   endif

   call modal_aero_wateruptake_init(shr_const_pi, errmsg_local, errflg_local)
   if (errflg_local /= 0) then
      call endrun('modal_aero_wateruptake_cam_init: ' // trim(errmsg_local))
   end if

   ! Register the diagnostic-list water uptake recompute with the aerosol
   ! interface (called by modal_aerosol_state%water_uptake for diagnostic
   ! radiation lists; wired at init because the portable modal aerosol
   ! schemes are not part of every build).
   call modal_aerosol_state_register_water_uptake_diag(modal_aero_wateruptake_diag)

end subroutine modal_aero_wateruptake_cam_init

!===============================================================================


subroutine modal_aero_wateruptake_dr(state, pbuf, aero_props, aero_state)
!-----------------------------------------------------------------------
!
! CAM specific driver for modal aerosol water uptake code.
!
!-----------------------------------------------------------------------

   use time_manager,  only: is_first_step
   use cam_history,   only: fieldname_len
   use tropopause,    only: tropopause_find_cam, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE

   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer
   class(aerosol_properties), intent(in), target :: aero_props
   class(aerosol_state), intent(in), target :: aero_state

   ! local variables

   integer  :: lchnk              ! chunk index
   integer  :: ncol               ! number of columns

   integer :: i, k, m
   integer :: itim_old
   integer :: nmodes
   integer :: tropLev(pcols)

   character(len=fieldname_len+3) :: fieldname

   real(r8), pointer :: h2ommr(:,:) ! specific humidity
   real(r8), pointer :: t(:,:)      ! temperatures (K)
   real(r8), pointer :: pmid(:,:)   ! layer pressure (Pa)

   real(r8), pointer :: cldn(:,:)      ! layer cloud fraction (0-1)
   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: dgncur_awet(:,:,:)
   real(r8), pointer :: wetdens(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)

   real(r8), pointer :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), pointer :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
   real(r8), pointer :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), pointer :: so4dryvol(:,:,:) ! single-particle-mean so4 dry volume (m3)
   real(r8), pointer :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
   real(r8), pointer :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)

   real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: specdens_1(:)

   real(r8), pointer :: sulfeq(:,:,:) ! H2SO4 equilibrium mixing ratios over particles (mol/mol)
   real(r8), allocatable :: sulfeq_local(:,:,:) ! local work array for _sub output
   real(r8), allocatable :: wtpct(:,:,:)  ! sulfate aerosol composition, weight % H2SO4
   real(r8), allocatable :: sulden(:,:,:) ! sulfate aerosol mass density (g/cm3)

   real(r8), allocatable :: alnsg(:)
   real(r8), pointer :: maer(:,:,:)      ! accumulated aerosol mode MRs

   real(r8) :: pm25(pcols,pver)      ! PM2.5 diagnostics
   real(r8) :: rhoair(pcols,pver)
   ! dmleung 20 Oct 2025 ++
   real(r8) :: pm25_mmr(pcols,pver)      ! PM2.5 mass mixing ratio     dmleung, 20 Nov 2023
   real(r8) :: pm1(pcols,pver)           ! PM1 mass conc
   real(r8) :: pm1_mmr(pcols,pver)       ! PM1 mass mixing ratio     dmleung, 20 Nov 2023
   real(r8) :: pm10(pcols,pver)          ! PM10 mass conc
   real(r8) :: pm10_mmr(pcols,pver)      ! PM10 mass mixing ratio     dmleung, 20 Nov 2023
   real(r8) :: pmtot_mmr(pcols,pver)      ! total PM mass mixing ratio
   ! dmleung --

   character(len=3) :: trnum       ! used to hold mode number (as characters)

   character(len=512) :: errmsg_local
   integer            :: errflg_local

   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   ! loop over all aerosol modes
   nmodes = aero_props%nbins()

   allocate( &
      wetrad(pcols,pver,nmodes),        &
      wetvol(pcols,pver,nmodes),        &
      wtrvol(pcols,pver,nmodes),        &
      wtpct(pcols,pver,nmodes),         &
      sulden(pcols,pver,nmodes),        &
      sulfeq_local(pcols,pver,nmodes),  &
      specdens_1(nmodes),              &
      alnsg(nmodes),                   &
      maer(pcols,pver,nmodes)          )

   call pbuf_get_field(pbuf, dgnum_idx,      dgncur_a )
   call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet )
   call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
   call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat)
   call pbuf_get_field(pbuf, hygro_idx,       hygro)
   call pbuf_get_field(pbuf, dryvol_idx,      dryvol)
   call pbuf_get_field(pbuf, dryrad_idx,      dryrad)
   call pbuf_get_field(pbuf, drymass_idx,     drymass)
   call pbuf_get_field(pbuf, so4dryvol_idx,   so4dryvol)
   call pbuf_get_field(pbuf, naer_idx,        naer)

   if (is_first_step()) then
      dgncur_awet(:,:,:) = dgncur_a(:,:,:)
   end if

   if (modal_strat_sulfate) then
      ! get tropopause level
      !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
      tropLev(:) = 0
      !REMOVECAM_END
      call tropopause_find_cam(state, tropLev, primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)
   endif

   h2ommr => state%q(:,:,1)
   t      => state%t
   pmid   => state%pmid

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   ! Zero output arrays (allocated at pcols, _sub only writes 1:ncol)
   wetrad(:,:,:)       = 0._r8
   wetvol(:,:,:)       = 0._r8
   wtrvol(:,:,:)       = 0._r8
   sulfeq_local(:,:,:) = 0._r8
   wtpct(:,:,:)        = 0._r8
   sulden(:,:,:)       = 0._r8
   maer(:,:,:)         = 0._r8

   call modal_aero_wateruptake_sub( &
      aero_props       = aero_props,            &
      aero_state       = aero_state,            &
      ncol             = ncol,                  &
      pver             = pver,                  &
      top_lev          = top_lev,               &
      do_strat_sulfate = modal_strat_sulfate,   &
      t                = t(:ncol,:),            &
      pmid             = pmid(:ncol,:),         &
      h2ommr           = h2ommr(:ncol,:),       &
      cldn             = cldn(:ncol,:),         &
      dryrad           = dryrad(:ncol,:,:),     &
      hygro            = hygro(:ncol,:,:),      &
      dryvol           = dryvol(:ncol,:,:),     &
      so4dryvol        = so4dryvol(:ncol,:,:),  &
      dgncur_awet      = dgncur_awet(:ncol,:,:), &
      troplev          = tropLev(:ncol),        &
      wetrad           = wetrad(:ncol,:,:),     &
      wetvol           = wetvol(:ncol,:,:),     &
      wtrvol           = wtrvol(:ncol,:,:),     &
      sulfeq           = sulfeq_local(:ncol,:,:), &
      wtpct            = wtpct(:ncol,:,:),      &
      sulden           = sulden(:ncol,:,:),     &
      specdens_1       = specdens_1,            &
      alnsg_out        = alnsg,                 &
      maer             = maer(:ncol,:,:),       &
      errmsg           = errmsg_local,          &
      errflg           = errflg_local)
   if (errflg_local /= 0) then
      call endrun('modal_aero_wateruptake_dr: ' // trim(errmsg_local))
   end if

   ! Copy sulfeq to pbuf and output strat sulfate diagnostics
   if (modal_strat_sulfate) then
      call pbuf_get_field(pbuf, sulfeq_idx, sulfeq)
      sulfeq(:,:,:) = sulfeq_local(:,:,:)

      do m = 1, nmodes
         fieldname = ' '
         write(fieldname,fmt='(a,i1)') 'wtpct_a',m
         call outfld(fieldname,wtpct(1:ncol,1:pver,m), ncol, lchnk )

         fieldname = ' '
         write(fieldname,fmt='(a,i1)') 'sulfeq_a',m
         call outfld(fieldname,sulfeq_local(1:ncol,1:pver,m), ncol, lchnk )

         fieldname = ' '
         write(fieldname,fmt='(a,i1)') 'sulden_a',m
         call outfld(fieldname,sulden(1:ncol,1:pver,m), ncol, lchnk )
      end do
   end if

   ! Post-processing: wet density, qaerwat, dgncur_awet update
   qaerwat = 0.0_r8

   do m = 1, nmodes

      do k = top_lev, pver
         do i = 1, ncol

            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
            qaerwat(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)

            ! compute aerosol wet density (kg/m3)
            if (wetvol(i,k,m) > 1.0e-30_r8) then
               wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m))/wetvol(i,k,m)
            else
               wetdens(i,k,m) = specdens_1(m)
            end if
         end do
      end do

   end do    ! modes

   ! Compute air density for PM diagnostics
   do k = top_lev, pver
      do i = 1, ncol
         rhoair(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

   pm25(:,:)=0._r8
   ! dmleung 20 Oct 2025 ++
   pm25_mmr(:,:)=0._r8
   pm1(:,:)=0._r8
   pm1_mmr(:,:)=0._r8
   pm10(:,:)=0._r8
   pm10_mmr(:,:)=0._r8
   pmtot_mmr(:,:)=0._r8
   ! dmleung --

   do m = 1, nmodes
      ! output to history
      write( trnum, '(i3.3)' ) m
      call outfld( 'wat_a'//trnum(3:3),  qaerwat(:,:,m),     pcols, lchnk)
      call outfld( 'dgnd_a'//trnum(2:3), dgncur_a(:,:,m),    pcols, lchnk)
      call outfld( 'dgnw_a'//trnum(2:3), dgncur_awet(:,:,m), pcols, lchnk)

      ! calculate PM2.5 diagnostics -- dgncur_a is zero above top_lev
      do k = top_lev, pver
         do i=1,ncol
            pm25(i,k) = pm25(i,k)+maer(i,k,m)*(1._r8-(0.5_r8 - 0.5_r8*erf(log(2.5e-6_r8/dgncur_a(i,k,m))/ &
                                              (2._r8**0.5_r8*alnsg(m)))))*rhoair(i,k)
            ! dmleung 20 Oct 2025: calculate other PM diagnostics ++
            pm25_mmr(i,k) = pm25_mmr(i,k)+maer(i,k,m)*(1._r8-(0.5_r8 - 0.5_r8*erf(log(2.5e-6_r8/dgncur_a(i,k,m))/ &
                                              (2._r8**0.5_r8*alnsg(m)))))    ! PM2.5 mass mixing ratio, dmleung
            pm1(i,k) = pm1(i,k)+maer(i,k,m)*(1._r8-(0.5_r8 - 0.5_r8*erf(log(1.0e-6_r8/dgncur_a(i,k,m))/ &
                                              (2._r8**0.5_r8*alnsg(m)))))*rhoair(i,k)
            pm1_mmr(i,k) = pm1_mmr(i,k)+maer(i,k,m)*(1._r8-(0.5_r8 - 0.5_r8*erf(log(1.0e-6_r8/dgncur_a(i,k,m))/ &
                                              (2._r8**0.5_r8*alnsg(m)))))    ! PM1 mass mixing ratio, dmleung
            pm10(i,k) = pm10(i,k)+maer(i,k,m)*(1._r8-(0.5_r8 - 0.5_r8*erf(log(10.0e-6_r8/dgncur_a(i,k,m))/ &
                                              (2._r8**0.5_r8*alnsg(m)))))*rhoair(i,k)
            pm10_mmr(i,k) = pm10_mmr(i,k)+maer(i,k,m)*(1._r8-(0.5_r8 - 0.5_r8*erf(log(10.0e-6_r8/dgncur_a(i,k,m))/ &
                                              (2._r8**0.5_r8*alnsg(m)))))    ! PM10 mass mixing ratio, dmleung
            pmtot_mmr(i,k) = pmtot_mmr(i,k)+maer(i,k,m)                        ! toal PM mass mixing ratio, dmleung
            ! dmleung --
         end do
      end do
   end do

   call outfld('PM25',     pm25(:,:),     pcols, lchnk)
   call outfld('PM25_SRF', pm25(:,pver),  pcols, lchnk)
   ! dmleung 20 Oct 2025 added history fields below ++
   call outfld('PM25_MMR', pm25_mmr(:,:), pcols, lchnk)
   call outfld('PM1_SRF',  pm1(:,pver), pcols, lchnk)
   call outfld('PM1_MMR',  pm1_mmr(:,:),  pcols, lchnk)
   call outfld('PM10_SRF', pm10(:,pver), pcols, lchnk)
   call outfld('PM10_MMR', pm10_mmr(:,:), pcols, lchnk)
   call outfld('PMTOT_MMR',pmtot_mmr(:,:),pcols, lchnk)
   call outfld('RHO_AIR',  rhoair(:,:),   pcols, lchnk)
   ! dmleung --

   deallocate(maer, alnsg)
   deallocate( &
      wetrad, wetvol, wtrvol, wtpct, sulden, sulfeq_local, specdens_1 )

end subroutine modal_aero_wateruptake_dr

!----------------------------------------------------------------------

   end module modal_aero_wateruptake_cam
