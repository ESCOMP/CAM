module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o, rair
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use wv_saturation,    only: qsat_water
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props
use cam_history,      only: addfld, add_default, outfld, horiz_only
use cam_logfile,      only: iulog
use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use phys_control,     only: phys_getopts
use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   modal_aero_wateruptake_init, &
   modal_aero_wateruptake_dr,   &
   modal_aero_wateruptake_sub,  &
   modal_aero_kohler

public :: modal_aero_wateruptake_reg

real(r8), parameter :: third = 1._r8/3._r8
real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8


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


logical, public :: modal_strat_sulfate = .false.   ! If .true. then MAM sulfate surface area density used in stratospheric heterogeneous chemistry

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_wateruptake_reg()

  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

   integer :: nmodes
   
   call rad_cnst_get_info(0, nmodes=nmodes)
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

subroutine modal_aero_wateruptake_init(pbuf2d)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field
   use infnan,       only : nan, assignment(=)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   real(r8) :: real_nan

   integer :: m, nmodes
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies

   character(len=3) :: trnum       ! used to hold mode number (as characters)
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

   call rad_cnst_get_info(0, nmodes=nmodes)

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
   
   call addfld('PM25',     (/ 'lev' /), 'A', 'kg/m3', 'PM2.5 concentration')
   call addfld('PM25_SRF', horiz_only,  'A', 'kg/m3', 'surface PM2.5 concentration')

   if (is_first_step()) then
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnumwet_idx, 0.0_r8)
      if (modal_strat_sulfate) then
      ! initialize fields in physics buffer to NaN (not a number) 
      ! so model will crash if used before initialization
         call pbuf_set_field(pbuf2d, sulfeq_idx, real_nan)
      endif
   endif

end subroutine modal_aero_wateruptake_init

!===============================================================================


subroutine modal_aero_wateruptake_dr(state, pbuf, list_idx_in, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m, hygro_m, dryvol_m, dryrad_m, drymass_m,&
                                     so4dryvol_m, naer_m)
!-----------------------------------------------------------------------
!
! CAM specific driver for modal aerosol water uptake code.
!
! *** N.B. *** The calculation has been enabled for diagnostic mode lists
!              via optional arguments.  If the list_idx arg is present then
!              all the optional args must be present.
!
!-----------------------------------------------------------------------

   use time_manager,  only: is_first_step
   use cam_history,   only: outfld, fieldname_len
   use tropopause,    only: tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer

   integer,  optional,          intent(in)    :: list_idx_in
   real(r8), optional,          pointer       :: dgnumdry_m(:,:,:)
   real(r8), optional,          pointer       :: dgnumwet_m(:,:,:)
   real(r8), optional,          pointer       :: qaerwat_m(:,:,:)
   real(r8), optional,          pointer       :: wetdens_m(:,:,:)
   real(r8), optional,          pointer       :: hygro_m(:,:,:)
   real(r8), optional,          pointer       :: dryvol_m(:,:,:)
   real(r8), optional,          pointer       :: dryrad_m(:,:,:)
   real(r8), optional,          pointer       :: drymass_m(:,:,:)
   real(r8), optional,          pointer       :: so4dryvol_m(:,:,:)
   real(r8), optional,          pointer       :: naer_m(:,:,:)

   ! local variables

   integer  :: lchnk              ! chunk index
   integer  :: ncol               ! number of columns
   integer  :: list_idx           ! radiative constituents list index
   integer  :: stat

   integer :: i, k, l, m
   integer :: itim_old
   integer :: nmodes
   integer :: nspec
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

   real(r8), pointer :: raer(:,:)        ! aerosol species MRs (kg/kg)
   real(r8), pointer :: maer(:,:,:)      ! accumulated aerosol mode MRs
   real(r8), pointer :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), pointer :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
   real(r8), pointer :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), pointer :: so4dryvol(:,:,:) ! single-particle-mean so4 dry volume (m3)
   real(r8), pointer :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
   real(r8), pointer :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)

   real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: rhcrystal(:)
   real(r8), allocatable :: rhdeliques(:)
   real(r8), allocatable :: specdens_1(:)

   real(r8), pointer :: sulfeq(:,:,:) ! H2SO4 equilibrium mixing ratios over particles (mol/mol)
   real(r8), allocatable :: wtpct(:,:,:)  ! sulfate aerosol composition, weight % H2SO4
   real(r8), allocatable :: sulden(:,:,:) ! sulfate aerosol mass density (g/cm3)
   
   real(r8) :: specdens, so4specdens
   real(r8) :: sigmag
   real(r8), allocatable :: alnsg(:)
   real(r8) :: rh(pcols,pver)        ! relative humidity (0-1)
   real(r8) :: dmean, qh2so4_equilib, wtpct_mode, sulden_mode

   real(r8) :: es(pcols)             ! saturation vapor pressure
   real(r8) :: qs(pcols)             ! saturation specific humidity

   real(r8) :: pm25(pcols,pver)      ! PM2.5 diagnostics     
   real(r8) :: rhoair(pcols,pver) 

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   character(len=32) :: spectype
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   list_idx = 0
   if (present(list_idx_in)) then
      list_idx = list_idx_in

      ! check that all optional args are present
      if (.not. present(dgnumdry_m) .or. .not. present(dgnumwet_m) .or. &
          .not. present(qaerwat_m)  .or. .not. present(wetdens_m)) then
         call endrun('modal_aero_wateruptake_dr called for'// &
                     'diagnostic list but required args not present')
      end if

      ! arrays for diagnostic calculations must be associated
      if (.not. associated(dgnumdry_m) .or. .not. associated(dgnumwet_m) .or. &
          .not. associated(qaerwat_m)  .or. .not. associated(wetdens_m)) then
         call endrun('modal_aero_wateruptake_dr called for'// &
                     'diagnostic list but required args not associated')
      end if

      if (modal_strat_sulfate) then
         call endrun('modal_aero_wateruptake_dr cannot be called with optional arguments and'// &
                     ' having modal_strat_sulfate set to true')
      end if
   end if

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (modal_strat_sulfate) then
     call pbuf_get_field(pbuf,  sulfeq_idx, sulfeq ) 
   endif

   allocate( &
      wetrad(pcols,pver,nmodes),   &
      wetvol(pcols,pver,nmodes),   &
      wtrvol(pcols,pver,nmodes),   &
      wtpct(pcols,pver,nmodes),    &
      sulden(pcols,pver,nmodes),   &
      rhcrystal(nmodes),           &
      rhdeliques(nmodes),          &
      specdens_1(nmodes),          &
      alnsg(nmodes)           )

   wtpct(:,:,:)     = 75._r8
   sulden(:,:,:)    = 1.923_r8

   if (list_idx == 0) then
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
   else
      dgncur_a    => dgnumdry_m
      dgncur_awet => dgnumwet_m
      qaerwat     => qaerwat_m
      wetdens     => wetdens_m
      hygro       => hygro_m
      dryvol      => dryvol_m
      dryrad      => dryrad_m
      drymass     => drymass_m
      so4dryvol   => so4dryvol_m
      naer        => naer_m
   end if

   if (modal_strat_sulfate) then
      ! get tropopause level
      call tropopause_find(state, tropLev, primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)
   endif

   h2ommr => state%q(:,:,1)
   t      => state%t
   pmid   => state%pmid

   allocate( maer(pcols,pver,nmodes))
   maer(:,:,:) = 0._r8

   do m = 1, nmodes

      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag, &
         rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      do l = 1, nspec

         ! accumulate the aerosol masses of each mode
         call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, raer)
         maer(:ncol,:,m)= maer(:ncol,:,m) + raer(:ncol,:)
            
         ! get species interstitial mixing ratio ('a')
         call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                     spectype=spectype)

         if (modal_strat_sulfate .and. (trim(spectype).eq.'sulfate')) then
            so4specdens=specdens
         end if

         if (l == 1) then
            ! save off these values to be used as defaults
            specdens_1(m) = specdens
         end if

      end do

      alnsg(m) = log(sigmag)

      if (modal_strat_sulfate) then
         do k = top_lev, pver
            do i = 1, ncol
               dmean = dgncur_awet(i,k,m)*exp(1.5_r8*alnsg(m)**2)
               call calc_h2so4_equilib_mixrat( t(i,k), pmid(i,k), h2ommr(i,k), dmean, &
                                               qh2so4_equilib, wtpct_mode, sulden_mode )
               sulfeq(i,k,m)  = qh2so4_equilib
               wtpct(i,k,m)   = wtpct_mode
               sulden(i,k,m)  = sulden_mode                        
            end do    ! i = 1, ncol
         end do    ! k = top_lev, pver
 
          fieldname = ' '
          write(fieldname,fmt='(a,i1)') 'wtpct_a',m
          call outfld(fieldname,wtpct(1:ncol,1:pver,m), ncol, lchnk )
 
          fieldname = ' '
          write(fieldname,fmt='(a,i1)') 'sulfeq_a',m
          call outfld(fieldname,sulfeq(1:ncol,1:pver,m), ncol, lchnk )
 
          fieldname = ' '
          write(fieldname,fmt='(a,i1)') 'sulden_a',m
          call outfld(fieldname,sulden(1:ncol,1:pver,m), ncol, lchnk )
 
       end if
      
    end do    ! m = 1, nmodes

   ! relative humidity calc

   itim_old    =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   do k = top_lev, pver
      call qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol))
      do i = 1, ncol
         if (qs(i) > h2ommr(i,k)) then
            rh(i,k) = h2ommr(i,k)/qs(i)
         else
            rh(i,k) = 0.98_r8
         endif
         rh(i,k) = max(rh(i,k), 0.0_r8)
         rh(i,k) = min(rh(i,k), 0.98_r8)
         if (cldn(i,k) .lt. 1.0_r8) then
            rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_r8 - cldn(i,k))  ! clear portion
         end if
         rh(i,k) = max(rh(i,k), 0.0_r8)
         rhoair(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rh, dryvol, so4dryvol, so4specdens, tropLev, &
      wetrad, wetvol, wtrvol, sulden, wtpct)

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

   if (list_idx == 0) then

      pm25(:,:)=0._r8

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
            end do
         end do
      end do

      call outfld('PM25',     pm25(:,:),    pcols, lchnk)
      call outfld('PM25_SRF', pm25(:,pver), pcols, lchnk)

   end if

   deallocate(maer, alnsg)
   deallocate( &
      wetrad, wetvol, wtrvol, wtpct, sulden, rhcrystal, rhdeliques, specdens_1 )
end subroutine modal_aero_wateruptake_dr

!===============================================================================

subroutine modal_aero_wateruptake_sub( &
   ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
   hygro, rh, dryvol, so4dryvol, so4specdens, troplev, &
   wetrad, wetvol, wtrvol, sulden, wtpct)

!-----------------------------------------------------------------------
!
! Purpose: Compute aerosol wet radius
!
! Method:  Kohler theory
!
! Author:  S. Ghan
!
!-----------------------------------------------------------------------

   
   ! Arguments
   integer, intent(in)  :: ncol                    ! number of columns
   integer, intent(in)  :: nmodes
   integer, intent(in)  :: troplev(:)

   real(r8), intent(in) :: rhcrystal(:)
   real(r8), intent(in) :: rhdeliques(:)
   real(r8), intent(in) :: dryrad(:,:,:)         ! dry volume mean radius of aerosol (m)
   real(r8), intent(in) :: hygro(:,:,:)          ! volume-weighted mean hygroscopicity (--)
   real(r8), intent(in) :: rh(:,:)               ! relative humidity (0-1)
   real(r8), intent(in) :: dryvol(:,:,:)         ! dry volume of single aerosol (m3)
   real(r8), intent(in) :: so4dryvol(:,:,:)      ! dry volume of sulfate in single aerosol (m3)
   real(r8), intent(in) :: so4specdens           ! mass density sulfate in single aerosol (kg/m3)
   real(r8), intent(in) :: wtpct(:,:,:)          ! sulfate aerosol composition, weight % H2SO4
   real(r8), intent(in) :: sulden(:,:,:)         ! sulfate aerosol mass density (g/cm3)

   real(r8), intent(out) :: wetrad(:,:,:)        ! wet radius of aerosol (m)
   real(r8), intent(out) :: wetvol(:,:,:)        ! single-particle-mean wet volume (m3)
   real(r8), intent(out) :: wtrvol(:,:,:)        ! single-particle-mean water volume in wet aerosol (m3)

   ! local variables

   integer :: i, k, m

   real(r8) :: hystfac                ! working variable for hysteresis
   !-----------------------------------------------------------------------


   ! loop over all aerosol modes
   do m = 1, nmodes

      hystfac = 1.0_r8 / max(1.0e-5_r8, (rhdeliques(m) - rhcrystal(m)))

      do k = top_lev, pver
         do i = 1, ncol

            if ( modal_strat_sulfate .and. (k<troplev(i))) then
               wetvol(i,k,m) = dryvol(i,k,m)-so4dryvol(i,k,m)
               wetvol(i,k,m) = wetvol(i,k,m)+so4dryvol(i,k,m)*so4specdens/sulden(i,k,m)/wtpct(i,k,m)/10._r8
               wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m))
               wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third
               wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m))
               wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m)
               wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_r8)
            else
              ! compute wet radius for each mode
              call modal_aero_kohler(dryrad(i:i,k,m), hygro(i:i,k,m), rh(i:i,k), wetrad(i:i,k,m), 1)

              wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m))
              wetvol(i,k,m) = pi43*wetrad(i,k,m)**3
              wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m))
              wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m)
              wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_r8)

              ! apply simple treatment of deliquesence/crystallization hysteresis
              ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
              ! the "upper curve" value, and the fraction is a linear function of rh
              if (rh(i,k) < rhcrystal(m)) then
                 wetrad(i,k,m) = dryrad(i,k,m)
                 wetvol(i,k,m) = dryvol(i,k,m)
                 wtrvol(i,k,m) = 0.0_r8
              else if (rh(i,k) < rhdeliques(m)) then
                 wtrvol(i,k,m) = wtrvol(i,k,m)*hystfac*(rh(i,k) - rhcrystal(m))
                 wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_r8)
                 wetvol(i,k,m) = dryvol(i,k,m) + wtrvol(i,k,m)
                 wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third
              end if
            end if

         end do  ! columns
      end do     ! levels

   end do ! modes

end subroutine modal_aero_wateruptake_sub

!-----------------------------------------------------------------------
      subroutine modal_aero_kohler(   &
          rdry_in, hygro, s, rwet_out, im )

! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35)

! for multiple aerosol types, assumes an internal mixture of aerosols

      implicit none

! arguments
      integer :: im         ! number of grid points to be processed
      real(r8) :: rdry_in(:)    ! aerosol dry radius (m)
      real(r8) :: hygro(:)      ! aerosol volume-mean hygroscopicity (--)
      real(r8) :: s(:)          ! relative humidity (1 = saturated)
      real(r8) :: rwet_out(:)   ! aerosol wet radius (m)

! local variables
      integer, parameter :: imax=200
      integer :: i, n, nsol

      real(r8) :: a, b
      real(r8) :: p40(imax),p41(imax),p42(imax),p43(imax) ! coefficients of polynomial
      real(r8) :: p30(imax),p31(imax),p32(imax) ! coefficients of polynomial
      real(r8) :: p
      real(r8) :: r3, r4
      real(r8) :: r(im)         ! wet radius (microns)
      real(r8) :: rdry(imax)    ! radius of dry particle (microns)
      real(r8) :: ss            ! relative humidity (1 = saturated)
      real(r8) :: slog(imax)    ! log relative humidity
      real(r8) :: vol(imax)     ! total volume of particle (microns**3)
      real(r8) :: xi, xr

      complex(r8) :: cx4(4,imax),cx3(3,imax)

      real(r8), parameter :: eps = 1.e-4_r8
      real(r8), parameter :: mw = 18._r8
      real(r8), parameter :: pi = 3.14159_r8
      real(r8), parameter :: rhow = 1._r8
      real(r8), parameter :: surften = 76._r8
      real(r8), parameter :: tair = 273._r8
      real(r8), parameter :: third = 1._r8/3._r8
      real(r8), parameter :: ugascon = 8.3e7_r8


!     effect of organics on surface tension is neglected
      a=2.e4_r8*mw*surften/(ugascon*tair*rhow)

      do i=1,im
           rdry(i) = rdry_in(i)*1.0e6_r8   ! convert (m) to (microns)
           vol(i) = rdry(i)**3          ! vol is r**3, not volume
           b = vol(i)*hygro(i)

!          quartic
           ss=min(s(i),1._r8-eps)
           ss=max(ss,1.e-10_r8)
           slog(i)=log(ss)
           p43(i)=-a/slog(i)
           p42(i)=0._r8
           p41(i)=b/slog(i)-vol(i)
           p40(i)=a*vol(i)/slog(i)
!          cubic for rh=1
           p32(i)=0._r8
           p31(i)=-b/a
           p30(i)=-vol(i)
      end do


       do 100 i=1,im

!       if(vol(i).le.1.e-20)then
        if(vol(i).le.1.e-12_r8)then
           r(i)=rdry(i)
           go to 100
        endif

        p=abs(p31(i))/(rdry(i)*rdry(i))
        if(p.lt.eps)then
!          approximate solution for small particles
           r(i)=rdry(i)*(1._r8+p*third/(1._r8-slog(i)*rdry(i)/a))
        else
           call makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1)
!          find smallest real(r8) solution
           r(i)=1000._r8*rdry(i)
           nsol=0
           do n=1,4
              xr=real(cx4(n,i))
              xi=aimag(cx4(n,i))
              if(abs(xi).gt.abs(xr)*eps) cycle  
              if(xr.gt.r(i)) cycle  
              if(xr.lt.rdry(i)*(1._r8-eps)) cycle  
              if(xr.ne.xr) cycle  
              r(i)=xr
              nsol=n
           end do  
           if(nsol.eq.0)then
              write(iulog,*)   &
               'ccm kohlerc - no real(r8) solution found (quartic)'
              write(iulog,*)'roots =', (cx4(n,i),n=1,4)
              write(iulog,*)'p0-p3 =', p40(i), p41(i), p42(i), p43(i)
              write(iulog,*)'rh=',s(i)
              write(iulog,*)'setting radius to dry radius=',rdry(i)
              r(i)=rdry(i)
!             stop
           endif
        endif

        if(s(i).gt.1._r8-eps)then
!          save quartic solution at s=1-eps
           r4=r(i)
!          cubic for rh=1
           p=abs(p31(i))/(rdry(i)*rdry(i))
           if(p.lt.eps)then
              r(i)=rdry(i)*(1._r8+p*third)
           else
              call makoh_cubic(cx3,p32,p31,p30,im)
!             find smallest real(r8) solution
              r(i)=1000._r8*rdry(i)
              nsol=0
              do n=1,3
                 xr=real(cx3(n,i))
                 xi=aimag(cx3(n,i))
                 if(abs(xi).gt.abs(xr)*eps) cycle  
                 if(xr.gt.r(i)) cycle  
                 if(xr.lt.rdry(i)*(1._r8-eps)) cycle  
                 if(xr.ne.xr) cycle  
                 r(i)=xr
                 nsol=n
              end do  
              if(nsol.eq.0)then
                 write(iulog,*)   &
                  'ccm kohlerc - no real(r8) solution found (cubic)'
                 write(iulog,*)'roots =', (cx3(n,i),n=1,3)
                 write(iulog,*)'p0-p2 =', p30(i), p31(i), p32(i)
                 write(iulog,*)'rh=',s(i)
                 write(iulog,*)'setting radius to dry radius=',rdry(i)
                 r(i)=rdry(i)
!                stop
              endif
           endif
           r3=r(i)
!          now interpolate between quartic, cubic solutions
           r(i)=(r4*(1._r8-s(i))+r3*(s(i)-1._r8+eps))/eps
        endif

  100 continue

! bound and convert from microns to m
      do i=1,im
         r(i) = min(r(i),30._r8) ! upper bound based on 1 day lifetime
         rwet_out(i) = r(i)*1.e-6_r8
      end do

      return
      end subroutine modal_aero_kohler


!-----------------------------------------------------------------------
      subroutine makoh_cubic( cx, p2, p1, p0, im )
!
!     solves  x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx)
      complex(r8) :: cx(3,imx)

      integer :: i
      real(r8) :: eps, q(imx), r(imx), sqrt3, third
      complex(r8) :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      save eps
      data eps/1.e-20_r8/

      third=1._r8/3._r8
      ci=cmplx(0._r8,1._r8,r8)
      sqrt3=sqrt(3._r8)
      cw=0.5_r8*(-1+ci*sqrt3)
      cwsq=0.5_r8*(-1-ci*sqrt3)

      do i=1,im
      if(p1(i).eq.0._r8)then
!        completely insoluble particle
         cx(1,i)=(-p0(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
      else
         q(i)=p1(i)/3._r8
         r(i)=p0(i)/2._r8
         crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
         crad(i)=sqrt(crad(i))

         cy(i)=r(i)-crad(i)
         if (abs(cy(i)).gt.eps) cy(i)=cy(i)**third
         cq=q(i)
         cz(i)=-cq/cy(i)

         cx(1,i)=-cy(i)-cz(i)
         cx(2,i)=-cw*cy(i)-cwsq*cz(i)
         cx(3,i)=-cwsq*cy(i)-cw*cz(i)
      endif
      enddo

      return
      end subroutine makoh_cubic


!-----------------------------------------------------------------------
      subroutine makoh_quartic( cx, p3, p2, p1, p0, im )

!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx), p3(imx)
      complex(r8) :: cx(4,imx)

      integer :: i
      real(r8) :: third, q(imx), r(imx)
      complex(r8) :: cb(imx), cb0(imx), cb1(imx),   &
                     crad(imx), cy(imx), czero


      czero=cmplx(0.0_r8,0.0_r8,r8)
      third=1._r8/3._r8

      do 10 i=1,im

      q(i)=-p2(i)*p2(i)/36._r8+(p3(i)*p1(i)-4*p0(i))/12._r8
      r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48._r8   &
       +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16

      crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
      crad(i)=sqrt(crad(i))

      cb(i)=r(i)-crad(i)
      if(cb(i).eq.czero)then
!        insoluble particle
         cx(1,i)=(-p1(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
         cx(4,i)=cx(1,i)
      else
         cb(i)=cb(i)**third

         cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6

         cb0(i)=sqrt(cy(i)*cy(i)-p0(i))
         cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i))

         cb(i)=p3(i)/2+cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i))
         crad(i)=sqrt(crad(i))
         cx(1,i)=(-cb(i)+crad(i))/2._r8
         cx(2,i)=(-cb(i)-crad(i))/2._r8

         cb(i)=p3(i)/2-cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i))
         crad(i)=sqrt(crad(i))
         cx(3,i)=(-cb(i)+crad(i))/2._r8
         cx(4,i)=(-cb(i)-crad(i))/2._r8
      endif
   10 continue

      return
      end subroutine makoh_quartic

!----------------------------------------------------------------------
      subroutine calc_h2so4_equilib_mixrat( temp, pres, qh2o, dmean, &
                                            qh2so4_equilib, wtpct, sulden )

      implicit none

      real(r8), intent(in)  :: temp            ! temperature (K)
      real(r8), intent(in)  :: pres            ! pressure (Pa)
      real(r8), intent(in)  :: qh2o            ! water vapor specific humidity (kg/kg)
      real(r8), intent(in)  :: dmean           ! mean diameter of particles in each mode
      real(r8), intent(out) :: qh2so4_equilib  ! h2so4 saturation mixing ratios over the particles (mol/mol)
      real(r8), intent(out) :: wtpct           ! sulfate composition, weight % H2SO4
      real(r8), intent(out) :: sulden          ! sulfate aerosol mass density (g/cm3)
      
      ! Local declarations
      real(r8)            :: qh2o_kelvin ! water vapor specific humidity adjusted for Kelvin effect (kg/kg)
      real(r8)            :: wtpct_flat  ! sulfate composition over a flat surface, weight % H2SO4
      real(r8)            :: fk1, fk4, fk4_1, fk4_2 
      real(r8)            :: factor_kulm                     ! Kulmala correction terms
      real(r8)            :: en, t, sig1, sig2, frac, surf_tens, surf_tens_mode, dsigma_dwt
      real(r8)            :: den1, den2, sulfate_density, drho_dwt
      real(r8)            :: akelvin, expon, akas, rkelvinH2O, rkelvinH2O_a, rkelvinH2O_b
      real(r8)            :: sulfequil, r
      real(r8), parameter :: t0_kulm     = 340._r8           !  T0 set in the low end of the Ayers measurement range (338-445K)
      real(r8), parameter :: t_crit_kulm = 905._r8           !  Critical temperature = 1.5 * Boiling point
      real(r8), parameter :: fk0         = -10156._r8 / t0_kulm + 16.259_r8   !  Log(Kulmala correction factor)
      real(r8), parameter :: fk2         = 1._r8 / t0_kulm
      real(r8), parameter :: fk3         = 0.38_r8 / (t_crit_kulm - t0_kulm)   
      real(r8), parameter :: RGAS = 8.31430e+07_r8           ! ideal gas constant (erg/mol/K)
      real(r8), parameter :: wtmol_h2so4 =  98.078479_r8     ! molecular weight of sulphuric acid
      real(r8), parameter :: wtmol_h2o   =  18.015280_r8     ! molecular weight of water vapor
      real(r8), parameter :: wtmol_air   =  28.9644_r8       ! molecular weight of dry air
      real(r8)            :: gv_wt_abcd_en(6,95), gvh2ovp(95)
      real(r8)            :: dnwtp(46), dnc0(46), dnc1(46)
      real(r8)            :: stwtp(15), stc0(15), stc1(15)
      integer             :: i, k

        
      data stwtp/0._r8, 23.8141_r8, 38.0279_r8, 40.6856_r8, 45.335_r8, 52.9305_r8, &
         56.2735_r8, 59.8557_r8, 66.2364_r8, 73.103_r8, 79.432_r8, 85.9195_r8, &
         91.7444_r8, 97.6687_r8, 100._r8/

      data stc0/117.564_r8, 103.303_r8, 101.796_r8, 100.42_r8, 98.4993_r8, 91.8866_r8, &
         88.3033_r8, 86.5546_r8, 84.471_r8, 81.2939_r8, 79.3556_r8, 75.608_r8, &
         70.0777_r8, 63.7412_r8, 61.4591_r8 /

      data stc1/-0.153641_r8, -0.0982007_r8, -0.0872379_r8, -0.0818509_r8,           &
         -0.0746702_r8, -0.0522399_r8, -0.0407773_r8, -0.0357946_r8, -0.0317062_r8,   &
         -0.025825_r8, -0.0267212_r8, -0.0269204_r8, -0.0276187_r8, -0.0302094_r8,    &
         -0.0303081_r8 /
           
       
      data dnwtp / 0._r8, 1._r8, 5._r8, 10._r8, 20._r8, 25._r8, 30._r8, 35._r8, 40._r8,  &
         41._r8, 45._r8, 50._r8, 53._r8, 55._r8, 56._r8, 60._r8, 65._r8, 66._r8, 70._r8, &
         72._r8, 73._r8, 74._r8, 75._r8, 76._r8, 78._r8, 79._r8, 80._r8, 81._r8, 82._r8, &
         83._r8, 84._r8, 85._r8, 86._r8, 87._r8, 88._r8, 89._r8, 90._r8, 91._r8, 92._r8, &
         93._r8, 94._r8, 95._r8, 96._r8, 97._r8, 98._r8, 100._r8 /
         
      data dnc0 / 1._r8, 1.13185_r8, 1.17171_r8, 1.22164_r8, 1.3219_r8, 1.37209_r8,         &
         1.42185_r8, 1.4705_r8, 1.51767_r8, 1.52731_r8, 1.56584_r8, 1.61834_r8, 1.65191_r8, &
         1.6752_r8, 1.68708_r8, 1.7356_r8, 1.7997_r8, 1.81271_r8, 1.86696_r8, 1.89491_r8,   &
         1.9092_r8, 1.92395_r8, 1.93904_r8, 1.95438_r8, 1.98574_r8, 2.00151_r8, 2.01703_r8, &
         2.03234_r8, 2.04716_r8, 2.06082_r8, 2.07363_r8, 2.08461_r8, 2.09386_r8, 2.10143_r8,&
         2.10764_r8, 2.11283_r8, 2.11671_r8, 2.11938_r8, 2.12125_r8, 2.1219_r8, 2.12723_r8, &
         2.12654_r8, 2.12621_r8, 2.12561_r8, 2.12494_r8, 2.12093_r8 /
         
      data dnc1 / 0._r8,  -0.000435022_r8, -0.000479481_r8, -0.000531558_r8, -0.000622448_r8, &
         -0.000660866_r8, -0.000693492_r8, -0.000718251_r8, -0.000732869_r8, -0.000735755_r8, &
         -0.000744294_r8, -0.000761493_r8, -0.000774238_r8, -0.00078392_r8, -0.000788939_r8,  &
         -0.00080946_r8, -0.000839848_r8, -0.000845825_r8, -0.000874337_r8, -0.000890074_r8,  &
         -0.00089873_r8, -0.000908778_r8, -0.000920012_r8, -0.000932184_r8, -0.000959514_r8,  &
         -0.000974043_r8, -0.000988264_r8, -0.00100258_r8, -0.00101634_r8, -0.00102762_r8,    &
         -0.00103757_r8, -0.00104337_r8, -0.00104563_r8, -0.00104458_r8, -0.00104144_r8,      &
         -0.00103719_r8, -0.00103089_r8, -0.00102262_r8, -0.00101355_r8, -0.00100249_r8,      &
         -0.00100934_r8, -0.000998299_r8, -0.000990961_r8, -0.000985845_r8, -0.000984529_r8,  &
         -0.000989315_r8 /  

      ! Saturation vapor pressure of sulfuric acid
      !  
      ! Limit extrapolation at extreme temperatures
      t=min(max(temp,140._r8),450._r8)
      
      !!  Calculate the weight % H2SO4 composition of sulfate
      call calc_h2so4_wtpct(t, pres, qh2o, wtpct_flat)
                  
      !!  Calculate surface tension (erg/cm2) of sulfate of 
      !!  different compositions as a linear function of temperature.
      i = 2 ! minimum wt% is 29.517
      do while (wtpct_flat.gt.stwtp(i))
       i = i + 1
      end do
      sig1 = stc0(i-1) + stc1(i-1) * t
      sig2 = stc0(i)   + stc1(i) * t
      ! calculate derivative needed later for kelvin factor for h2o      
      dsigma_dwt = (sig2-sig1) / (stwtp(i)-stwtp(i-1))
      surf_tens = sig1 + dsigma_dwt*(wtpct_flat-stwtp(i))    
      
      !!  Calculate density (g/cm3) of sulfate of 
      !!  different compositions as a linear function of temperature.
      i = 6 ! minimum wt% is 29.517
      do while (wtpct_flat .gt. dnwtp(i))
        i = i + 1
      end do
      den1 = dnc0(i-1) + dnc1(i-1) * t
      den2 = dnc0(i)   + dnc1(i) * t
      ! Calculate derivative needed later for Kelvin factor for H2O      
      drho_dwt = (den2-den1) / (dnwtp(i)-dnwtp(i-1))
      sulfate_density = den1 + drho_dwt * (wtpct_flat-dnwtp(i-1))

      r=dmean*100._r8/2._r8 ! calcuate mode radius (cm) from diameter (m)

      ! Adjust for Kelvin effect for water
      rkelvinH2O_b = 1 + wtpct_flat * drho_dwt / &
         sulfate_density - 3._r8 * wtpct_flat * dsigma_dwt / (2._r8*surf_tens)

      rkelvinH2O_a = 2._r8 * wtmol_h2so4 * surf_tens / &
           (sulfate_density * RGAS * t * r)     

      rkelvinH2O = exp (rkelvinH2O_a*rkelvinH2O_b)

      qh2o_kelvin = qh2o/rkelvinH2O
      call calc_h2so4_wtpct(t, pres, qh2o_kelvin, wtpct)


      wtpct=max(wtpct,wtpct_flat)

      ! Parameterized fit to Giauque's (1959) enthalpies v. wt %:
      en = 4.184_r8 * (23624.8_r8 - 1.14208e8_r8 / ((wtpct - 105.318_r8)**2 + 4798.69_r8))
      en = max(en, 0.0_r8)

      !!  Calculate surface tension (erg/cm2) of sulfate of 
      !!  different compositions as a linear function of temperature.
      i = 2 ! minimum wt% is 29.517
      do while (wtpct.gt.stwtp(i))
       i=i+1
      end do
      sig1=stc0(i-1)+stc1(i-1)*t
      sig2=stc0(i)+stc1(i)*t
      frac=(stwtp(i)-wtpct)/(stwtp(i)-stwtp(i-1))
      surf_tens_mode=sig1*frac+sig2*(1.0_r8-frac)     

      !!  Calculate density (g/cm3) of sulfate of 
      !!  different compositions as a linear function of temperature.
      i = 6 ! minimum wt% is 29.517
      do while (wtpct .gt. dnwtp(i))
        i=i+1
      end do
      den1=dnc0(i-1)+dnc1(i-1)*t
      den2=dnc0(i)+dnc1(i)*t
      frac=(dnwtp(i)-wtpct)/(dnwtp(i)-dnwtp(i-1))
      sulden=den1*frac+den2*(1.0_r8-frac)

      ! Ayers' (1980) fit to sulfuric acid equilibrium vapor pressure:
      ! (Remember this is the log)
      ! SULFEQ=-10156/Temp+16.259-En/(8.3143*Temp)
      !
      ! Kulmala correction (J. CHEM. PHYS. V.93, No.1, 1 July 1990)
      fk1   = -1._r8 / t
      fk4_1 = log(t0_kulm / t)
      fk4_2 = t0_kulm / t
      fk4   = 1.0_r8 + fk4_1 - fk4_2
      factor_kulm = fk1 + fk2 + fk3 * fk4

      ! This is for pure H2SO4
      sulfequil = fk0 + 10156._r8 * factor_kulm

      ! Adjust for WTPCT composition:
      sulfequil = sulfequil - en / (8.3143_r8 * t)

      ! Take the exponential:
      sulfequil = exp(sulfequil)

      ! Convert atmospheres ==> Pa
      sulfequil = sulfequil * 1.01325e5_r8  

      ! Convert Pa ==> mol/mol
      sulfequil = sulfequil / pres

      ! Calculate Kelvin curvature factor for H2SO4 interactively with temperature:
      ! (g/mol)*(erg/cm2)/(K * g/cm3 * erg/mol/K) = cm
      akelvin = 2._r8 * wtmol_h2so4 * surf_tens_mode / (t * sulden * RGAS)

      expon = akelvin / r  ! divide by mode radius (cm) 
      expon = max(-100._r8, expon)
      expon = min(100._r8, expon)
      akas = exp( expon )
      qh2so4_equilib = sulfequil * akas ! reduce H2SO4 equilibrium mixing ratio by Kelvin curvature factor

      return
      end subroutine calc_h2so4_equilib_mixrat


!----------------------------------------------------------------------
      subroutine calc_h2so4_wtpct( temp, pres, qh2o, wtpct )
      
  !!  This function calculates the weight % H2SO4 composition of 
  !!  sulfate aerosol, using Tabazadeh et. al. (GRL, 1931, 1997).
  !!  Rated for T=185-260K, activity=0.01-1.0
  !!
  !!  Argument list input:   
  !!    temp = temperature (K)
  !!    pres = atmospheric pressure (Pa)
  !!    qh2o = water specific humidity (kg/kg)
  !!
  !!  Output:
  !!    wtpct = weight % H2SO4 in H2O/H2SO4 particle (0-100)
  !!
  !! @author Mike Mills
  !! @ version October 2013

      use wv_saturation, only: qsat_water
      
      implicit none

      real(r8), intent(in)  :: temp  ! temperature (K)
      real(r8), intent(in)  :: pres  ! pressure (Pa)
      real(r8), intent(in)  :: qh2o  ! water vapor specific humidity (kg/kg)
      real(r8), intent(out) :: wtpct ! sulfate weight % H2SO4 composition
      
      ! Local declarations
      real(r8)            :: atab1,btab1,ctab1,dtab1,atab2,btab2,ctab2,dtab2
      real(r8)            :: contl, conth, contt, conwtp
      real(r8)            :: activ
      real(r8)            :: es ! saturation vapor pressure over water (Pa) (dummy)
      real(r8)            :: qs ! saturation specific humidity over water (kg/kg)

      ! calculate saturation specific humidity over pure water, qs (kg/kg)
      call qsat_water(temp, pres, es, qs)
      
      !  Activity = water specific humidity (kg/kg) / equilibrium water (kg/kg)
      activ = qh2o/qs
       
      if (activ.lt.0.05_r8) then
        activ = max(activ,1.e-6_r8)    ! restrict minimum activity
        atab1 	= 12.37208932_r8	
        btab1 	= -0.16125516114_r8
        ctab1 	= -30.490657554_r8
        dtab1 	= -2.1133114241_r8
        atab2 	= 13.455394705_r8	
        btab2 	= -0.1921312255_r8
        ctab2 	= -34.285174607_r8
        dtab2 	= -1.7620073078_r8
      elseif (activ.ge.0.05_r8.and.activ.le.0.85_r8) then
        atab1 	= 11.820654354_r8
        btab1 	= -0.20786404244_r8
        ctab1 	= -4.807306373_r8
        dtab1 	= -5.1727540348_r8
        atab2 	= 12.891938068_r8	
        btab2 	= -0.23233847708_r8
        ctab2 	= -6.4261237757_r8
        dtab2 	= -4.9005471319_r8
      elseif (activ.gt.0.85_r8) then
        activ = min(activ,1._r8)      ! restrict maximum activity
        atab1 	= -180.06541028_r8
        btab1 	= -0.38601102592_r8
        ctab1 	= -93.317846778_r8
        dtab1 	= 273.88132245_r8
        atab2 	= -176.95814097_r8
        btab2 	= -0.36257048154_r8
        ctab2 	= -90.469744201_r8
        dtab2 	= 267.45509988_r8
      else
        write(iulog,*) 'calc_h2so4_wtpct: invalid activity: activ,qh2o,qs,temp,pres=',activ,qh2o,qs,temp,pres
        call endrun( 'calc_h2so4_wtpct error' )
        return
      endif

      contl = atab1*(activ**btab1)+ctab1*activ+dtab1
      conth = atab2*(activ**btab2)+ctab2*activ+dtab2

      contt = contl + (conth-contl) * ((temp -190._r8)/70._r8)
      conwtp = (contt*98._r8) + 1000._r8

      wtpct = (100._r8*contt*98._r8)/conwtp
      wtpct = min(max(wtpct,25._r8),100._r8) ! restrict between 1 and 100 %
      
      return
      end subroutine calc_h2so4_wtpct


!----------------------------------------------------------------------

   end module modal_aero_wateruptake


