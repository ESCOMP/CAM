module crmclouds_camaerosols
#if (defined CRM)
#if (defined MODAL_AERO) 
!---------------------------------------------------------------------------------------------
! Purpose: 
! 
!  Provides the necessary subroutines to use cloud fields from the CRM model to drive the 
!  aerosol-related subroutines in CAM. Several taskes:
!     i) to fill the physics buffers with those diagnosed from the CRM clouds.  
!    ii) to provide the interface for some physics prcoesses, such as droplet activaiton, 
!         and convetive transport. 
!
!  An alternative (and better?) approach is to use the ECPP (explicit-cloud parameterized-pollutant). 
!  This will be done later.
!
!  Revision history: 
!  July, 27, 2009: Minghuai Wang
! 
!-------------------------------------------------------------------------------------------- 
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use cam_abortutils,  only: endrun
   use crmdims,         only: crm_nx, crm_ny, crm_nz
   use physics_buffer,  only: physics_buffer_desc, pbuf_get_field, pbuf_get_index, pbuf_old_tim_idx
   use physics_types,   only: physics_state, physics_state_copy, physics_ptend
   use ref_pres,        only: top_lev => clim_modal_aero_top_lev
   use wv_saturation,   only: qsat_water
   implicit none
   private
   save

   public :: spcam_modal_aero_wateruptake_dr
   public :: crmclouds_mixnuc_tend 
   public :: crmclouds_diag 
   public :: crmclouds_convect_tend

!======================================================================================================
contains 

subroutine spcam_modal_aero_wateruptake_dr(state,pbuf)

!-----------------------------------------------------------------------
!
! SPCAM specific driver for modal aerosol water uptake code.
!
!-----------------------------------------------------------------------

   use time_manager,          only: is_first_step
   use modal_aero_wateruptake,only: modal_aero_wateruptake_sub
   use physconst,             only: pi, rhoh2o
   use rad_constituents,      only: rad_cnst_get_info, rad_cnst_get_mode_props, rad_cnst_get_aer_props


   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer

   ! local variables

   real(r8), parameter :: third = 1._r8/3._r8
   real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8

   integer  :: ncol               ! number of columns

   integer :: i, k, m
   integer :: nmodes
   integer :: nspec
   integer :: mm

   integer :: dgnumwet_idx, qaerwat_idx, wetdens_ap_idx, cld_idx

   integer :: dgnum_idx      = 0
   integer :: hygro_idx      = 0
   integer :: dryvol_idx     = 0
   integer :: dryrad_idx     = 0
   integer :: drymass_idx    = 0
   integer :: so4dryvol_idx  = 0
   integer :: naer_idx       = 0

   real(r8), allocatable :: wtrvol_grid(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)
   real(r8), allocatable :: wetvol_grid(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: ncount_clear(:,:,:)   ! to count the fraction of clear sky part

   real(r8), pointer :: h2ommr_crm(:,:,:,:)  ! specfic humidity in CRM domain
   real(r8), pointer :: t_crm(:,:,:,:)  ! temperature at the CRM domain
   real(r8), pointer :: cldn_crm(:,:,:,:)  ! cloud fraction in CRM domain
   real(r8), pointer :: qaerwat_crm(:, :, :, :, :)  ! aerosol water at CRM domain
   real(r8), pointer :: dgncur_awet_crm(:, :, :, :, :)   ! wet mode diameter at CRM domain

   real(r8),allocatable :: es_crm(:)         ! saturation vapor pressure
   real(r8),allocatable :: qs_crm(:)         ! saturation specific humidity
   real(r8),allocatable :: cldnt(:,:)        ! temporal variables
   real(r8),allocatable :: rh_crm(:,:,:,:)   ! Relative humidity at the CRM grid
   real(r8),allocatable :: specdens_1(:)

   real(r8),pointer :: dgncur_a(:,:,:)
   real(r8),pointer :: drymass(:,:,:)
   real(r8),pointer :: dryrad(:,:,:)


   real(r8), pointer :: dgncur_awet(:,:,:)
   real(r8), pointer :: wetdens(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)

   real(r8), pointer :: h2ommr(:,:) ! specific humidity
   real(r8), pointer :: t(:,:)      ! temperatures (K)
   real(r8), pointer :: pmid(:,:)   ! layer pressure (Pa)
   real(r8), pointer :: cldn(:,:)   ! layer cloud fraction (0-1)

   real(r8), allocatable :: maer(:,:,:)      ! aerosol wet mass MR (including water) (kg/kg-air)
   real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)
   real(r8), allocatable :: wtpct(:,:,:)     ! sulfate aerosol composition, weight % H2SO4
   real(r8), allocatable :: sulden(:,:,:)    ! sulfate aerosol mass density (g/cm3)

   real(r8), pointer :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), pointer :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
   real(r8), pointer :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), pointer :: so4dryvol(:,:,:) ! dry volume of sulfate in single aerosol (m3)

   real(r8) :: specdens, so4specdens
   integer               :: troplev(pcols)

   real(r8), allocatable :: rhcrystal(:)
   real(r8), allocatable :: rhdeliques(:)

   real(r8) :: es(pcols)             ! saturation vapor pressure
   real(r8) :: qs(pcols)             ! saturation specific humidity



   real(r8) :: rh(pcols,pver)        ! relative humidity (0-1)


   real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)

   integer :: ii, jj, l
   integer :: idx
   integer :: itim_old


   !-----------------------------------------------------------------------

   ncol = state%ncol

   call rad_cnst_get_info(0, nmodes=nmodes)

   allocate(&
      es_crm(pcols),                           &
      qs_crm(pcols),                           &
      cldnt(pcols, pver),                      &
      rh_crm(pcols, crm_nx, crm_ny, pver), & 
      wtrvol_grid(pcols,pver,nmodes),          &
      wetvol_grid(pcols,pver,nmodes),          &
      ncount_clear(pcols,pver,nmodes),         &
      dgncur_a(pcols,pver,nmodes),             &
      drymass(pcols,pver,nmodes),              &
      specdens_1(nmodes)                       )

   allocate(  &
      wetrad(pcols,pver,nmodes),   & 
      wetvol(pcols,pver,nmodes),   &
      wtrvol(pcols,pver,nmodes),   &
      wtpct(pcols,pver,nmodes),    &
      sulden(pcols,pver,nmodes),   &
      rhcrystal(nmodes),           &
      rhdeliques(nmodes)           )

   wtpct(:,:,:)     = 75._r8
   sulden(:,:,:)    = 1.923_r8

   dgnum_idx       = pbuf_get_index('DGNUM')
   hygro_idx       = pbuf_get_index('HYGRO')
   dryvol_idx      = pbuf_get_index('DRYVOL')
   dryrad_idx      = pbuf_get_index('DRYRAD')
   drymass_idx     = pbuf_get_index('DRYMASS')
   so4dryvol_idx   = pbuf_get_index('SO4DRYVOL')
   naer_idx        = pbuf_get_index('NAER')
   dgnumwet_idx    = pbuf_get_index('DGNUMWET')
   qaerwat_idx     = pbuf_get_index('QAERWAT')
   wetdens_ap_idx  = pbuf_get_index('WETDENS_AP')
   cld_idx         = pbuf_get_index('CLD')


   idx = pbuf_get_index('CRM_QV_RAD')
   call pbuf_get_field (pbuf, idx, h2ommr_crm)
   idx = pbuf_get_index('CRM_T_RAD')
   call pbuf_get_field (pbuf, idx, t_crm)
   idx = pbuf_get_index('CRM_CLD_RAD')
   call pbuf_get_field (pbuf, idx, cldn_crm)
   idx = pbuf_get_index('CRM_QAERWAT')
   call pbuf_get_field (pbuf, idx, qaerwat_crm)
   idx = pbuf_get_index('CRM_DGNUMWET')
   call pbuf_get_field (pbuf, idx, dgncur_awet_crm)

   ncount_clear  = 0.0_r8
   wtrvol_grid   = 0.0_r8
   wetvol_grid   = 0.0_r8

   call pbuf_get_field(pbuf, hygro_idx,       hygro)
   call pbuf_get_field(pbuf, dryvol_idx,      dryvol)
   call pbuf_get_field(pbuf, dryrad_idx,      dryrad)
   call pbuf_get_field(pbuf, drymass_idx,     drymass)
   call pbuf_get_field(pbuf, so4dryvol_idx,   so4dryvol)
   call pbuf_get_field(pbuf, naer_idx,        naer)

   call pbuf_get_field(pbuf, dgnum_idx,       dgncur_a )
   call pbuf_get_field(pbuf, dgnumwet_idx,    dgncur_awet )
   call pbuf_get_field(pbuf, wetdens_ap_idx,  wetdens)
   call pbuf_get_field(pbuf, qaerwat_idx,     qaerwat)

   dgncur_awet(:,:,:) = dgncur_a(:,:,:)
   qaerwat            = 0._r8 

   h2ommr => state%q(:,:,1)
   t      => state%t
   pmid   => state%pmid

   do m = 1, nmodes
      ! get mode properties
      call rad_cnst_get_mode_props(0, m, rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))
      ! get mode info
      call rad_cnst_get_info(0, m, nspec=nspec)

      do l = 1, nspec

         ! get species interstitial mixing ratio ('a')
         call rad_cnst_get_aer_props(0, m, l, density_aer=specdens)

         if (l == 1) then
            ! save off these values to be used as defaults
            specdens_1(m) = specdens
         end if

      end do

   end do

   itim_old    =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   do jj = 1, crm_ny
      do ii = 1, crm_nx
         do k = top_lev, pver
            mm=pver-k+1
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
            end do

            if (mm <= crm_nz) call qsat_water(t_crm(:ncol,ii,jj,mm), &
                                                 pmid(:ncol,k), es_crm(:ncol), qs_crm(:ncol))
            do i = 1, ncol
               rh_crm(i, ii, jj, k) = rh(i,k)
               if(mm.le.crm_nz) then
                  rh_crm(i, ii, jj, k) = h2ommr_crm(i,ii,jj,mm)/qs_crm(i)
                  rh_crm(i, ii, jj, k) = max(rh_crm(i, ii, jj, k), 0.0_r8)
                  rh_crm(i, ii, jj, k) = min(rh_crm(i, ii, jj, k), 0.98_r8)
                  if(cldn_crm(i, ii, jj, mm).gt.0.5_r8) then
                     ! aerosol water uptake is not calculaed at overcast sky in MMF
                     rh_crm(i, ii, jj, k) = 0.0_r8
                  end if
                end if

                rh(i,k) = rh_crm(i, ii, jj, k)
                cldnt(i, k) = cldn(i,k)
                mm=pver-k+1
                if(mm.le.crm_nz) then
                   cldnt(i,k) = cldn_crm(i, ii, jj, mm)
                end if

                do m=1,nmodes
                   ncount_clear(i,k,m) = ncount_clear(i,k,m) + (1._r8 - cldnt(i,k))
                end do
            end do
         end do

         call modal_aero_wateruptake_sub( &
                                    ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
                                    hygro, rh, dryvol, so4dryvol, so4specdens, tropLev, &
                                    wetrad, wetvol, wtrvol, sulden, wtpct)
         do m = 1, nmodes
            do k = top_lev, pver
               do i = 1, ncol
                  dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
                  if(k.ge.pver-crm_nz+1) then
                     qaerwat_crm(i,ii,jj,pver-k+1,m) = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)
                     dgncur_awet_crm(i,ii,jj,pver-k+1,m) = dgncur_awet(i,k,m)
                  end if
                  wtrvol_grid(i,k,m) = wtrvol_grid(i,k,m) + wtrvol(i,k,m)*(1._r8-cldnt(i,k))
                  wetvol_grid(i,k,m) = wetvol_grid(i,k,m) + wetvol(i,k,m)*(1._r8-cldnt(i,k))
                  qaerwat(i,k,m) = qaerwat(i,k,m)+ rhoh2o*naer(i,k,m)*wtrvol(i,k,m) * (1-cldnt(i,k))

               end do
            end do
         end do
      end do
   end do

   do m = 1, nmodes
      do k = 1, pver
         do i = 1, ncol

            if(ncount_clear(i,k,m).gt.1.0e-10_r8) then
              qaerwat(i,k,m) = qaerwat(i,k,m)/ncount_clear(i,k,m)
              wetvol_grid(i,k,m)=wetvol_grid(i,k,m)/ncount_clear(i,k,m)
              wtrvol_grid(i,k,m)=wtrvol_grid(i,k,m)/ncount_clear(i,k,m)
              if (wetvol_grid(i,k,m) > 1.0e-30_r8) then
                 wetdens(i,k,m) = (drymass(i,k,m) + &
                                   rhoh2o*wtrvol_grid(i,k,m))/wetvol_grid(i,k,m)
              else
                 wetdens(i,k,m) = specdens_1(m)
              end if
              wetrad(i,k,m) = max(dryrad(i,k,m), (wetvol_grid(i,k,m)/pi43)**third)
              dgncur_awet(i,k,m) = dgncur_a(i,k,m)*   &
                 (wetrad(i,k,m)/dryrad(i,k,m))
            else
              dgncur_awet(i,k,m) = dgncur_a(i,k,m)
              qaerwat(i,k,m) = 0.0_r8
              wetdens(i,k,m) = specdens_1(m)
            end if
         end do   ! ncol
      end do   ! pver
   end do   ! nmodes



   deallocate(&
      es_crm,           &
      qs_crm,           &
      cldnt,            &
      rh_crm,           &
      wtrvol_grid,      &
      wetvol_grid,      &
      ncount_clear      )

   deallocate(wetrad, wetvol, wtrvol, wtpct, sulden, rhcrystal, rhdeliques, specdens_1)

end subroutine spcam_modal_aero_wateruptake_dr


!------------------------------------------------------------------------------------------------------
subroutine crmclouds_mixnuc_tend (state, ptend, dtime, cflx, pblht, pbuf,   &
                   wwqui_cen, wwqui_cloudy_cen, wwqui_bnd, wwqui_cloudy_bnd )
!-----------------------------------------------------------------------------------------------------
!
! Purpose: to calculate aerosol tendency from dropelt activation and mixing. 
!          Adopted from mmicro_pcond in cldwat2m.F90
!
!------------------------------------------------------------------------------------------------------
  use physics_types,    only: physics_state, physics_ptend, physics_tend, physics_ptend_init
  use physics_buffer,   only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
  use physconst,        only: gravit, rair, karman
  use constituents,     only: cnst_get_ind, pcnst, cnst_species_class, cnst_spec_class_gas
  use time_manager,     only: is_first_step
  use cam_history,      only: outfld
  use ndrop,         only: dropmixnuc
  use modal_aero_data
  use rad_constituents, only: rad_cnst_get_info

! Input 
  type(physics_state), intent(in)    :: state   ! state variables
  type(physics_buffer_desc), pointer :: pbuf(:)
  real(r8), intent(in) :: pblht(pcols)          ! PBL height (meter)
  real(r8), intent(in)  :: dtime                ! timestep
  real(r8), intent(in) :: cflx(pcols,pcnst)     ! constituent flux from surface
  real(r8), intent(in) :: wwqui_cen(pcols, pver)           ! vertical velocity variance in quiescent class (m2/s2)
  real(r8), intent(in) :: wwqui_cloudy_cen(pcols, pver)    ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
  real(r8), intent(in) :: wwqui_bnd(pcols, pver+1)         ! vertical velocity variance in quiescent class (m2/s2)
  real(r8), intent(in) :: wwqui_cloudy_bnd(pcols, pver+1)  ! vertical velocity variance in quiescent, and cloudy class (m2/s2)

! output
  type(physics_ptend), intent(out) :: ptend   ! package tendencies

! Local variables
  integer i,k,m, k1, k2
  integer ifld, itim
  integer ixcldliq, ixcldice, ixnumliq
  integer l,lnum,lnumcw,lmass,lmasscw
  integer :: lchnk                  ! chunk identifier
  integer :: ncol                   ! number of atmospheric columns
  integer :: nmodes
 
  real(r8) :: nc(pcols, pver)       ! droplet number concentration (#/kg)
  real(r8) :: nctend(pcols, pver)   ! change in droplet number concentration
  real(r8) :: omega(pcols, pver)    ! grid-averaaged vertical velocity 
  real(r8) :: qc(pcols, pver)       ! liquid water content (kg/kg)
  real(r8) :: qi(pcols, pver)       ! ice water content (kg/kg) 
  real(r8) :: lcldn(pcols, pver)
  real(r8) :: lcldo(pcols, pver) 
  real(r8) :: cldliqf(pcols, pver) 

  real(r8) :: wsub(pcols, pver)     ! subgrid vertical velocity
  real(r8) :: ekd_crm(pcols, pverp)  ! diffusivity
  real(r8) :: kkvh_crm(pcols, pverp)  ! eddy diffusivity
  real(r8) :: zs(pcols, pver)       ! inverse of distance between levels (meter)
  real(r8) :: dz(pcols, pver)       ! layer depth (m)
  real(r8) :: cs(pcols, pver)       ! air density
  real(r8) :: lc(pcols, pverp)       ! mixing length (m)
  real(r8) :: zheight(pcols, pverp)   ! height at lay interface (m)
  
  real(r8) :: alc(pcols, pverp)        ! asymptotic length scale (m)
  real(r8) :: tendnd(pcols, pver)      ! tendency of cloud droplet number concentrations (not used in the MMF) 
 
  real(r8),allocatable :: factnum(:,:,:)  ! activation fraction for aerosol number

  real(r8) :: qcld, qsmall

  logical :: dommf=.true.              ! value insignificant, if present, means that dropmixnuc is called the mmf part. 

! Variables in the physics buffer:
  real(r8), pointer, dimension(:,:) :: cldn    ! cloud fractin at the current time step
  real(r8), pointer, dimension(:,:) :: cldo   ! cloud fraction at the previous time step
  real(r8), pointer, dimension(:,:) :: acldy_cen ! liquid cloud fraction at the previous time step from ECPP
  real(r8), pointer, dimension(:,:) ::  kkvh    ! vertical diffusivity
  real(r8), pointer, dimension(:,:) :: tke          ! turbulence kenetic energy 
  real(r8), pointer, dimension(:,:) :: tk_crm     ! m2/s

  logical :: lq(pcnst)

  lchnk = state%lchnk
  ncol  = state%ncol

  qsmall = 1.e-18_r8

  call rad_cnst_get_info(0, nmodes=nmodes)
  allocate(factnum(pcols,pver,nmodes))

  lq(:) = .false.
  do m=1,ntot_amode
    lnum=numptr_amode(m)
    if(lnum>0)then
       lq(lnum)= .true.
    endif
    do l=1,nspec_amode(m)
      lmass=lmassptr_amode(l,m)
      lq(lmass)= .true.
    enddo
  enddo

  call physics_ptend_init(ptend,state%psetcols,'crmclouds_mixnuc', lq=lq)

!
! In the MMF model, turbulent mixing for tracer species are turned off in tphysac.
! So the turbulent for gas species mixing are added here.
!
   do m=1, pcnst
      if(cnst_species_class(m).eq.cnst_spec_class_gas) then
        ptend%lq(m) = .true.
      end if
   end do

  itim = pbuf_old_tim_idx ()
  ifld = pbuf_get_index ('CLD')
  call pbuf_get_field(pbuf, ifld, cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  ifld = pbuf_get_index ('CLDO')
  call pbuf_get_field(pbuf, ifld, cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  ifld = pbuf_get_index ('ACLDY_CEN')
  call pbuf_get_field(pbuf, ifld, acldy_cen)
  ifld = pbuf_get_index('kvh')
  call pbuf_get_field(pbuf, ifld, kkvh)

  ifld=pbuf_get_index('tke')
  call pbuf_get_field(pbuf, ifld, tke)

  ifld = pbuf_get_index('TK_CRM')
  call pbuf_get_field(pbuf, ifld, tk_crm)


  if (is_first_step()) then
     kkvh(:,:)= 0.0_r8
     tke(:,:) = 0.0_r8
  endif

  do i=1, ncol
    do k=1, pver-1
    zs(i,k) = 1._r8/(state%zm(i,k)-state%zm(i,k+1))
    end do
    zs(i,pver) = zs(i,pver-1)

! calculate height at layer interface (simple calculation)
    zheight(i,pverp) = 0.0_r8
    do k=pver, 1, -1
      zheight(i,k) = zheight(i,k+1) + state%pdel(i,k)/state%pmid(i,k)*(rair*state%t(i,k)/gravit)
    end do

! calculate mixing length
! from Holtslag and Boville, 1993, J. Climate. 
!
    do k=1, pverp
      if(zheight(i,k).le.pblht(i)) then
        alc(i,k) = 300._r8
      else
        alc(i,k) = 30._r8+270._r8*exp(1._r8-zheight(i,k)/pblht(i))
      endif
      lc(i,k) = alc(i,k)*karman*zheight(i,k)/(alc(i,k)+karman*zheight(i,k))
    enddo 
  end do

  call outfld('LENGC', lc, pcols, lchnk)

  kkvh_crm = 0._r8
  do i=1, ncol
    do k=1, pver

! from vertical variance in the quiescent class, which excldues 
! the contribution from strong updraft and downdraft. 
       wsub(i,k) = sqrt(wwqui_cloudy_cen(i,k))  ! use variance in cloudy quiescent area
       wsub(i,k) = min(wsub(i,k), 10._r8) 
       wsub(i,k) = max(0.20_r8, wsub(i,k))
    end do   ! end k

    do k=1, pver+1

      k1=min(k, pver)
      k2=max(k-1, 1)
!
! calculate ekd_crm from wsub in the cloudy quiescent class (following a part of ndrop.F90)
      ekd_crm(i,k) = min(10.0_r8, max(0.20_r8, sqrt(wwqui_cloudy_bnd(i,k))))* lc(i,k) 
      kkvh_crm(i,k) = ekd_crm(i,k)

! set kkvh to kkvh_crm so this will be used in dropmixnuc in the mmf call
      kkvh(i,k) = kkvh_crm(i,k)

    end do   !end k

  end do

  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)
  call cnst_get_ind('NUMLIQ', ixnumliq)

  qc(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
  qi(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
  nc(:ncol,:pver) = state%q(:ncol,:pver,ixnumliq)
  cldliqf(:,:)    = 1._r8
  lcldn(:,:)      = 0._r8
  lcldo(:,:)      = 0._r8
  

  do k=1,pver
   do i=1,ncol
      qcld=qc(i,k)+qi(i,k)
      if(qcld.gt.qsmall)then

#ifdef ECPP
!
!      When ECPP is called, activation associated with cloud fraction change is treated in ECPP.
!      so set two cloud fractio be the same here. 
!      But ECPP still did not treat activation associated with turbulent scale motion, and is
!      done in dropmixnuc
         lcldn(i,k)=acldy_cen(i,k)
         lcldo(i,k)=acldy_cen(i,k)
#else
         lcldn(i,k)=cldn(i,k)*qc(i,k)/qcld
         lcldo(i,k)=cldo(i,k)*qc(i,k)/qcld
#endif
      else
         lcldn(i,k)=0._r8
         lcldo(i,k)=0._r8
      endif
    enddo
  enddo

! should we set omega to be zero ??
  omega(:ncol, :) = state%omega(:ncol, :)

  call dropmixnuc(state, ptend, dtime, pbuf, wsub, lcldn, lcldo, cldliqf, tendnd, factnum, dommf )

! this part is moved into tphysbc after aerosol stuffs. 
!

  deallocate(factnum)

end subroutine crmclouds_mixnuc_tend
!======================================================================================================

!------------------------------------------------------------------------------------------------------
subroutine crmclouds_convect_tend(state,  ptend,  ztodt,  pbuf)
!-----------------------------------------------------------------
!
! Purpose: to do convective transport of tracer species using the cloud fields from CRM and using the 
!          subroutine of convtran. 
!
! Minghuai Wang, July, 2009: adopted from zm_conv_tend_2 
!
!------------------------------------------------------------------------------------------------------
   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
   use constituents,  only: pcnst, cnst_get_ind
   use zm_conv,       only: convtran
   use error_messages, only: alloc_err  

! Arguments
! Input variables:
   type(physics_state), intent(in ) :: state          ! Physics state variables
   real(r8), intent(in) :: ztodt

! Output variables:
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer

! Local variables
   integer :: i, lchnk, istat
   integer :: ncol
   integer :: nstep
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
   real(r8), dimension(pcols,pver) :: dpdry
   real(r8), dimension(pcols,pver) :: dp   ! layer thickness in mbs (between upper/lower interface).
   real(r8), dimension(pcols) :: dsubcld   !  wg layer thickness in mbs between lcl and maxi.

! physics buffer fields 
   integer itim, ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble


   real(r8), pointer, dimension(:,:) :: mu  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: eu  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: du  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: md  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: ed  !(pcols,pver,begchunk:endchunk)

   real(r8), pointer, dimension(:) :: jtr8   !(pcols,begchunk:endchunk)
        ! wg top  level index of deep cumulus convection.
   real(r8), pointer, dimension(:) :: maxgr8 !(pcols,begchunk:endchunk)
        ! wg gathered values of maxi.
   real(r8), pointer, dimension(:) :: ideepr8 !(pcols,begchunk:endchunk)               
        ! w holds position of gathered points vs longitude index

   integer :: jt(pcols)
   integer :: maxg(pcols)
   integer :: ideep(pcols) 
   integer :: lengath !(begchunk:endchunk)
   logical :: lq(pcnst)

!
! Initialize
!

   lq(:) = .true.
   lq(1)        = .false.
   lq(ixcldice) = .false.
   lq(ixcldliq) = .false.

   call physics_ptend_init(ptend,state%psetcols,'convtran2',lq=lq)

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols,pver,pcnst/) )

   ifld = pbuf_get_index('MU_CRM')
   call pbuf_get_field(pbuf, ifld, mu)
   ifld = pbuf_get_index('MD_CRM')
   call pbuf_get_field(pbuf, ifld, md)
   ifld = pbuf_get_index('DU_CRM')
   call pbuf_get_field(pbuf, ifld, du)
   ifld = pbuf_get_index('EU_CRM')
   call pbuf_get_field(pbuf, ifld, eu)
   ifld = pbuf_get_index('ED_CRM')
   call pbuf_get_field(pbuf, ifld, ed)
   ifld = pbuf_get_index('JT_CRM')
   call pbuf_get_field(pbuf, ifld, jtr8)
   ifld = pbuf_get_index('MX_CRM')
   call pbuf_get_field(pbuf, ifld, maxgr8)
   ifld = pbuf_get_index('IDEEP_CRM')
   call pbuf_get_field(pbuf, ifld, ideepr8)


! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk
  ncol = state%ncol

   nstep = get_nstep()

!
!     Convective transport of all trace species except cloud liquid 
!     and cloud ice done here because we need to do the scavenging first
!     to determine the interstitial fraction.
!
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   
! Is this ok to get the index???
   jt = int(jtr8+0.5_r8)
   maxg = int(maxgr8+0.5_r8)
   ideep = int(ideepr8+0.5_r8)

! calculate lengath from ideep
   lengath = 0
   do i=1, ncol
    if(ideep(i).ge.1) then
      lengath = lengath + 1
    endif
   end do

!
! initialize dpdry for call to convtran 
! it is used for tracers of dry smixing ratio type
!
   dpdry = 0._r8
   do i = 1,lengath
     dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
     dp(i,:) = state%pdel(ideep(i),:)/100._r8
   end do

! dsubdld is not used in convtran, and is set to be zero. 
   dsubcld = 0._r8


   call convtran (lchnk,                                        &
                  ptend%lq,state%q, pcnst,  mu(:,:), md(:,:),   &
                  du(:,:), eu(:,:), ed(:,:), dp(:,:), dsubcld(:),  &
                  jt(:),maxg(:),ideep(:), 1, lengath,  &
                  nstep,   fracis,  ptend%q, dpdry, ztodt  )

end subroutine crmclouds_convect_tend
!=====================================================================================================

!------------------------------------------------------------------------------------------------------
subroutine crmclouds_diag

end subroutine crmclouds_diag
!======================================================================================================

#endif
#endif /*CRM*/

end module crmclouds_camaerosols
