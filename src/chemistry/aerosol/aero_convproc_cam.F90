module aero_convproc_cam
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to aerosol/trace-gas convective cloud processing scheme
!
! currently these routines assume stratiform and convective clouds only interact
! through the detrainment of convective cloudborne material into stratiform clouds
!
! thus the stratiform-cloudborne aerosols (in the qqcw array) are not processed
! by the convective up/downdrafts, but are affected by the detrainment
!
! Author: R. C. Easter
!
!---------------------------------------------------------------------------------

   use shr_kind_mod,    only: r8=>shr_kind_r8
   use shr_kind_mod,    only: shr_kind_cs

   use spmd_utils,      only: masterproc
   use physconst,       only: gravit, rair
   use physconst,       only: pi, rhoh2o, rh2o, latvap, cpair
   use ppgrid,          only: pver, pcols
   use constituents,    only: pcnst, cnst_get_ind
   use constituents,    only: cnst_species_class, cnst_spec_class_aerosol
   use phys_control,    only: phys_getopts

   use physics_types,   only: physics_state, physics_ptend
   use physics_buffer,  only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
   use cam_history,     only: outfld, addfld, add_default, horiz_only
   use cam_logfile,     only: iulog
   use cam_abortutils,  only: endrun

   use aerosol_properties_mod, only: aerosol_properties
   use aerosol_state_mod, only: aerosol_state, ptr2d_t

   use aero_convproc,   only: aero_convproc_run
   use aero_convproc,   only: use_cwaer_for_activate_maxsat, convproc_method_activate
   use aero_convproc,   only: method1_activate_nlayers, method2_activate_smaxmax
   use aero_convproc,   only: method_reduce_actfrac, factor_reduce_actfrac

   implicit none
   private

   public :: aero_convproc_readnl
   public :: aero_convproc_init
   public :: aero_convproc_intr

! namelist options
! NOTE: These are the defaults for CAM6.
   logical, protected, public :: deepconv_wetdep_history = .true.
   logical, protected, public :: convproc_do_deep = .true.
! NOTE: These are the defaults for the Eaton/Wang parameterization.
   logical, protected, public :: convproc_do_evaprain_atonce = .false.
   real(r8), protected, public    :: convproc_pom_spechygro = -1._r8
   real(r8), protected, public    :: convproc_wup_max       = 4.0_r8

   logical, parameter :: apply_convproc_tend_to_ptend = .true.

   logical :: convproc_do_aer

! physics buffer indices
   integer :: fracis_idx         = 0

   integer :: rprddp_idx         = 0
   integer :: nevapr_dpcu_idx    = 0

   integer :: icwmrdp_idx        = 0
   integer :: dp_frac_idx        = 0

   integer :: zm_eu_idx          = 0
   integer :: zm_du_idx          = 0
   integer :: zm_ed_idx          = 0
   integer :: zm_dp_idx          = 0
   integer :: zm_jt_idx          = 0
   integer :: zm_maxg_idx        = 0
   integer :: zm_ideep_idx       = 0

   integer :: istat

   integer :: nbins = 0
   integer :: ncnstaer = 0

   integer, allocatable :: aer_cnst_ndx(:)

   character(len=32), allocatable :: cnst_name_extd(:,:) ! (2,ncnstaer)

contains

!=========================================================================================
   subroutine aero_convproc_readnl(nlfile)

      use namelist_utils, only: find_group_name
      use spmd_utils,     only: mpicom, masterprocid, mpi_real8, mpi_logical

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'aero_convproc_readnl'

      namelist /aerosol_convproc_opts/ deepconv_wetdep_history, convproc_do_deep, &
         convproc_do_evaprain_atonce, convproc_pom_spechygro, convproc_wup_max

      ! Read namelist
      if (masterproc) then
         open( newunit=unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'aerosol_convproc_opts', status=ierr)
         if (ierr == 0) then
            read(unitn, aerosol_convproc_opts, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
      end if

      ! Broadcast namelist variables
      call mpi_bcast( deepconv_wetdep_history,  1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast( convproc_do_deep,  1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast( convproc_do_evaprain_atonce,  1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast( convproc_pom_spechygro,  1, mpi_real8, masterprocid, mpicom, ierr)
      call mpi_bcast( convproc_wup_max,  1, mpi_real8, masterprocid, mpicom, ierr)

      if (masterproc) then
         write(iulog,*) subname//': deepconv_wetdep_history = ',deepconv_wetdep_history
         write(iulog,*) subname//': convproc_do_deep = ',convproc_do_deep
         write(iulog,*) subname//': convproc_do_evaprain_atonce = ',convproc_do_evaprain_atonce
         write(iulog,*) subname//': convproc_pom_spechygro = ',convproc_pom_spechygro
         write(iulog,*) subname//': convproc_wup_max = ', convproc_wup_max
      end if

   end subroutine aero_convproc_readnl

!=========================================================================================

   subroutine aero_convproc_init(aero_props)

      class(aerosol_properties), intent(in) :: aero_props

      integer :: m, mm, l, ndx, astat
      integer :: npass_calc_updraft
      logical :: history_aerosol
      character(len=32) :: name_a, name_c

      character(len=*), parameter :: prefix = 'aero_convproc_init: '

      nbins = aero_props%nbins()
      ncnstaer = aero_props%ncnst_tot()

      allocate(aer_cnst_ndx(ncnstaer),stat=astat)
      if (astat/=0) then
         call endrun(prefix//'aer_cnst_ndx allocation error')
      end if
      allocate(cnst_name_extd(2,ncnstaer),stat=astat)
      if (astat/=0) then
         call endrun(prefix//'cnst_name_extd allocation error')
      end if

      aer_cnst_ndx(:) = -1

      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            if (l==0) then
               call aero_props%num_names(m, name_a, name_c)
            else
               call aero_props%mmr_names(m,l, name_a, name_c)
            endif
            cnst_name_extd(1,mm) = name_a
            cnst_name_extd(2,mm) = name_c

            call cnst_get_ind(trim(name_a), ndx, abort=.false.)
            aer_cnst_ndx(mm) = ndx
         end do
      end do

      call phys_getopts( history_aerosol_out=history_aerosol, &
                        convproc_do_aer_out = convproc_do_aer )

      call addfld('DP_MFUP_MAX', horiz_only, 'A', 'kg/m2', &
                  'Deep conv. column-max updraft mass flux' )
      call addfld('DP_WCLDBASE', horiz_only, 'A', 'm/s', &
                  'Deep conv. cloudbase vertical velocity' )
      call addfld('DP_KCLDBASE', horiz_only, 'A', '1', &
                  'Deep conv. cloudbase level index' )

      ! output wet deposition fields to history
      !    I = in-cloud removal;     E = precip-evap resuspension
      !    C = convective (total);   D = deep convective
      ! note that the precip-evap resuspension includes that resulting from
      !    below-cloud removal, calculated in mz_aero_wet_intr
      if (convproc_do_aer .and. apply_convproc_tend_to_ptend ) then

         do m = 1, aero_props%nbins()
            do l = 0, aero_props%nmasses(m)
               mm = aero_props%indexer(m,l)

               ndx = aer_cnst_ndx(mm)

               if ( deepconv_wetdep_history ) then
                  call addfld (trim(cnst_name_extd(1,mm))//'SFSID', &
                               horiz_only,  'A','kg/m2/s','Wet deposition flux (incloud, deep convective) at surface')
                  call addfld (trim(cnst_name_extd(1,mm))//'SFSED', &
                               horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, deep convective) at surface')
                  if (history_aerosol) then
                     call add_default(trim(cnst_name_extd(1,mm))//'SFSID', 1, ' ')
                     call add_default(trim(cnst_name_extd(1,mm))//'SFSED', 1, ' ')
                  end if
               end if

            end do
         end do
      end if

      if ( history_aerosol .and. convproc_do_aer ) then
         call add_default( 'DP_MFUP_MAX', 1, ' ' )
         call add_default( 'DP_WCLDBASE', 1, ' ' )
         call add_default( 'DP_KCLDBASE', 1, ' ' )
      end if

      fracis_idx      = pbuf_get_index('FRACIS')

      rprddp_idx      = pbuf_get_index('RPRDDP')
      nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')

      icwmrdp_idx     = pbuf_get_index('ICWMRDP')
      dp_frac_idx     = pbuf_get_index('DP_FRAC')

      zm_eu_idx       = pbuf_get_index('ZM_EU')
      zm_du_idx       = pbuf_get_index('ZM_DU')
      zm_ed_idx       = pbuf_get_index('ZM_ED')
      zm_dp_idx       = pbuf_get_index('ZM_DP')
      zm_jt_idx       = pbuf_get_index('ZM_JT')
      zm_maxg_idx     = pbuf_get_index('ZM_MAXG')
      zm_ideep_idx    = pbuf_get_index('ZM_IDEEP')

      if (masterproc ) then

         write(iulog,'(a,l12)')     'aero_convproc_init - convproc_do_aer               = ', &
            convproc_do_aer
         write(iulog,'(a,l12)')     'aero_convproc_init - use_cwaer_for_activate_maxsat = ', &
            use_cwaer_for_activate_maxsat
         write(iulog,'(a,l12)')     'aero_convproc_init - apply_convproc_tend_to_ptend  = ', &
            apply_convproc_tend_to_ptend
         write(iulog,'(a,i12)')     'aero_convproc_init - convproc_method_activate      = ', &
            convproc_method_activate
         write(iulog,'(a,i12)')     'aero_convproc_init - method1_activate_nlayers      = ', &
            method1_activate_nlayers
         write(iulog,'(a,1pe12.4)') 'aero_convproc_init - method2_activate_smaxmax      = ', &
            method2_activate_smaxmax
         write(iulog,'(a,i12)')     'aero_convproc_init - method_reduce_actfrac         = ', &
            method_reduce_actfrac
         write(iulog,'(a,1pe12.4)') 'aero_convproc_init - factor_reduce_actfrac         = ', &
            factor_reduce_actfrac

         npass_calc_updraft = 1
         if ( (method_reduce_actfrac == 2)      .and. &
             (factor_reduce_actfrac >= 0.0_r8) .and. &
             (factor_reduce_actfrac <= 1.0_r8) ) npass_calc_updraft = 2
         write(iulog,'(a,i12)')     'aero_convproc_init - npass_calc_updraft            = ', &
            npass_calc_updraft

      end if

   end subroutine aero_convproc_init

!=========================================================================================

   subroutine aero_convproc_intr( aero_props, aero_state, state, ptend, pbuf, ztodt,             &
                                 nsrflx_mzaer2cnvpr, qsrflx_mzaer2cnvpr,  &
                                 aerdepwetis, dcondt_resusp3d )
!-----------------------------------------------------------------------
!
! Convective cloud processing (transport, activation/resuspension,
!    wet removal) of aerosols and trace gases.
!    (Currently no aqueous chemistry and no trace-gas wet removal)
! Does aerosols    when convproc_do_aer is .true.
!
! Does deep convection
! Uses mass fluxes, cloud water, precip production from the
!    convective cloud routines
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

      ! Arguments
      class(aerosol_properties), intent(in) :: aero_props
      class(aerosol_state), intent(in) :: aero_state

      type(physics_state),target,intent(in )   :: state      ! Physics state variables
      type(physics_ptend),       intent(inout) :: ptend      ! %lq set in aero_model_wetdep
      type(physics_buffer_desc), pointer       :: pbuf(:)
      real(r8), intent(in) :: ztodt                          ! model physics timestep [s]

      integer,  intent(in)    :: nsrflx_mzaer2cnvpr
      real(r8), intent(in)    :: qsrflx_mzaer2cnvpr(pcols,ncnstaer,nsrflx_mzaer2cnvpr)
      real(r8), intent(inout) :: aerdepwetis(pcols,pcnst)  ! aerosol wet deposition (interstitial)
      real(r8), intent(inout) :: dcondt_resusp3d(ncnstaer,pcols,pver)

      ! Local variables
      integer, parameter :: nsrflx = 5        ! last dimension of qsrflx
      integer  :: l, m, mm, ndx, lchnk
      integer  :: ncol

      real(r8) :: dqdt(pcols,pver,ncnstaer)
      real(r8) :: dt

      real(r8) :: q(pcols,pver,ncnstaer)
      real(r8) :: qsrflx(pcols,ncnstaer,nsrflx)
      real(r8), pointer :: qptr(:,:)

      real(r8) :: sflxic(pcols,ncnstaer)
      real(r8) :: sflxid(pcols,ncnstaer)
      real(r8) :: sflxec(pcols,ncnstaer)
      real(r8) :: sflxed(pcols,ncnstaer)

      type(ptr2d_t) :: raer(ncnstaer)     ! aerosol mass, number mixing ratios
      type(ptr2d_t) :: qqcw(ncnstaer)

      logical  :: applytend

      !-------------------------------------------------------------------------------------------------

      ! Initialize
      lchnk = state%lchnk
      ncol  = state%ncol
      dt = ztodt

      sflxic(:,:) = 0.0_r8
      sflxid(:,:) = 0.0_r8
      sflxec(:,:) = 0.0_r8
      sflxed(:,:) = 0.0_r8

      call aero_state%get_states( aero_props, raer, qqcw )

      ! prepare for deep conv processing
      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)

            mm = aero_props%indexer(m,l)
            ndx = aer_cnst_ndx(mm)

            sflxec(1:ncol,mm) = qsrflx_mzaer2cnvpr(1:ncol,mm,1)
            sflxed(1:ncol,mm) = qsrflx_mzaer2cnvpr(1:ncol,mm,2)

            applytend = .false.
            if ( ndx > 0 ) then
               applytend = ptend%lq(ndx)
            endif

            qptr => raer(mm)%fld

            if ( applytend ) then
               ! calc new q (after calcaersize and mz_aero_wet_intr)
               q(1:ncol,:,mm) = max( 0.0_r8, qptr(1:ncol,:) + dt*ptend%q(1:ncol,:,ndx) )
            else
               ! use old q
               q(1:ncol,:,mm) = qptr(1:ncol,:)
            end if

         end do
      end do

      dqdt(:,:,:) = 0.0_r8
      qsrflx(:,:,:) = 0.0_r8

      if (convproc_do_aer) then

         ! do deep conv processing
         if (convproc_do_deep) then
            call aero_convproc_dp_intr( aero_props, &
                                       state, pbuf, dt,                          &
                                       q, dqdt, nsrflx, qsrflx, dcondt_resusp3d )

            ! apply deep conv processing tendency

            do m = 1, aero_props%nbins()
               do l = 0, aero_props%nmasses(m)
                  mm = aero_props%indexer(m,l)
                  ndx = aer_cnst_ndx(mm)

                  if ( apply_convproc_tend_to_ptend ) then
                     ! add dqdt onto ptend%q and set ptend%lq
                     if (ndx>0) then ! advected species
                        ptend%q(1:ncol,:,ndx) = ptend%q(1:ncol,:,ndx) + dqdt(1:ncol,:,mm)
                     else
                        raer(mm)%fld(1:ncol,:) = max( 0.0_r8, raer(mm)%fld(1:ncol,:) + dqdt(1:ncol,:,mm) * dt )
                     end if
                  end if

                  ! these used for history file wetdep diagnostics
                  sflxic(1:ncol,mm) = sflxic(1:ncol,mm) + qsrflx(1:ncol,mm,4)
                  sflxid(1:ncol,mm) = sflxid(1:ncol,mm) + qsrflx(1:ncol,mm,4)
                  sflxec(1:ncol,mm) = sflxec(1:ncol,mm) + qsrflx(1:ncol,mm,5)
                  sflxed(1:ncol,mm) = sflxed(1:ncol,mm) + qsrflx(1:ncol,mm,5)

                  ! this used for surface coupling
                  if (ndx>0) then
                     aerdepwetis(1:ncol,ndx) = aerdepwetis(1:ncol,ndx) &
                                               + qsrflx(1:ncol,mm,4) + qsrflx(1:ncol,mm,5)
                  end if
               end do
            end do

         end if

      end if ! (convproc_do_aer) then

      if (convproc_do_aer .and. apply_convproc_tend_to_ptend ) then

         do m = 1, aero_props%nbins()
            do l = 0, aero_props%nmasses(m)
               mm = aero_props%indexer(m,l)

               ndx = aer_cnst_ndx(mm)

               if (ndx>0) call outfld( trim(cnst_name_extd(1,mm))//'SFWETC', aerdepwetis(:,ndx), pcols, lchnk )
               call outfld( trim(cnst_name_extd(1,mm))//'SFSIC', sflxic(:,mm), pcols, lchnk )
               call outfld( trim(cnst_name_extd(1,mm))//'SFSEC', sflxec(:,mm), pcols, lchnk )

               if ( deepconv_wetdep_history ) then
                  call outfld( trim(cnst_name_extd(1,mm))//'SFSID', sflxid(:,mm), pcols, lchnk )
                  call outfld( trim(cnst_name_extd(1,mm))//'SFSED', sflxed(:,mm), pcols, lchnk )
               end if
            end do
         end do

      end if

   end subroutine aero_convproc_intr

!=========================================================================================

   subroutine aero_convproc_dp_intr( aero_props,  &
                                    state, pbuf, dt,                          &
                                    q, dqdt, nsrflx, qsrflx,  dcondt_resusp3d)
!-----------------------------------------------------------------------
!
! Convective cloud processing (transport, activation/resuspension,
!    wet removal) of aerosols and trace gases.
!    (Currently no aqueous chemistry and no trace-gas wet removal)
! Does aerosols    when convproc_do_aer is .true.
!
! This routine does deep convection
! Uses mass fluxes, cloud water, precip production from the
!    convective cloud routines
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

      ! Arguments
      class(aerosol_properties), intent(in) :: aero_props

      type(physics_state),       intent(in ) :: state          ! Physics state variables
      type(physics_buffer_desc), pointer     :: pbuf(:)

      real(r8), intent(in) :: dt                         ! delta t (model time increment)

      real(r8), intent(in)    :: q(pcols,pver,ncnstaer)
      real(r8), intent(inout) :: dqdt(pcols,pver,ncnstaer)
      integer,  intent(in)    :: nsrflx
      real(r8), intent(inout) :: qsrflx(pcols,ncnstaer,nsrflx)
      real(r8), intent(inout) :: dcondt_resusp3d(ncnstaer,pcols,pver)

      integer :: i, l, m, mm
      integer :: lchnk

      real(r8) :: dpdry(pcols,pver)     ! layer delta-p-dry (mb)
      real(r8) :: fracice(pcols,pver)   ! Ice fraction of cloud droplets
      real(r8) :: xx_mfup_max(pcols), xx_wcldbase(pcols), xx_kcldbase(pcols)

      ! updraft interface TMR + wet-deposition TMR tendency diagnostics returned
      !    from aero_convproc_run for the WETC/CONU history fields
      real(r8) :: conu2(pcols,pver,2,ncnstaer)
      real(r8) :: dcondt2(pcols,pver,2,ncnstaer)

      character(len=512) :: errmsg
      integer            :: errflg

      ! physics buffer fields
      real(r8), pointer :: fracis(:,:,:)  ! fraction of transported species that are insoluble
      real(r8), pointer :: rprddp(:,:)    ! Deep conv precip production (kg/kg/s - grid avg)
      real(r8), pointer :: evapcdp(:,:)   ! Deep conv precip evaporation (kg/kg/s - grid avg)
      real(r8), pointer :: icwmrdp(:,:)   ! Deep conv cloud condensate (kg/kg - in cloud)
      real(r8), pointer :: dp_frac(:,:)   ! Deep conv cloud frac (0-1)

      ! deep conv variables
      real(r8), pointer :: du(:,:)        ! Mass detrain rate from updraft (pcols,pver)
      real(r8), pointer :: eu(:,:)        ! Mass entrain rate into updraft (pcols,pver)
      real(r8), pointer :: ed(:,:)        ! Mass entrain rate into downdraft (pcols,pver)
      ! eu, ed, du are "d(massflux)/dp" and are all positive
      real(r8), pointer :: dp(:,:)        ! Delta pressure between interfaces (pcols,pver)
      integer,  pointer :: jt(:)          ! Index of cloud top for each column (pcols)
      integer,  pointer :: maxg(:)        ! Index of cloud bottom for each column (pcols)
      integer,  pointer :: ideep(:)       ! Gathering array (pcols)
      integer           :: lengath        ! Gathered min lon indices over which to operate

      ! Initialize

      lchnk = state%lchnk

      ! Associate pointers with physics buffer fields
      call pbuf_get_field(pbuf, rprddp_idx,      rprddp)
      call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp)
      call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp)
      call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac)
      call pbuf_get_field(pbuf, fracis_idx,      fracis)
      call pbuf_get_field(pbuf, zm_eu_idx,       eu)
      call pbuf_get_field(pbuf, zm_du_idx,       du)
      call pbuf_get_field(pbuf, zm_ed_idx,       ed)
      call pbuf_get_field(pbuf, zm_dp_idx,       dp)
      call pbuf_get_field(pbuf, zm_jt_idx,       jt)
      call pbuf_get_field(pbuf, zm_maxg_idx,     maxg)
      call pbuf_get_field(pbuf, zm_ideep_idx,    ideep)

      lengath = count(ideep > 0)

      fracice(:,:) = 0.0_r8

      ! initialize dpdry (units=mb), which is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1, lengath
         dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
      end do

      call aero_convproc_run( aero_props, 'deep', lchnk,   dt,      &
                             state%t,    state%pmid, q, du,      eu,      &
                             ed,         dp,         dpdry,      jt,      &
                             maxg,       ideep,      1,          lengath, &
                             dp_frac,    icwmrdp,    rprddp,     evapcdp, &
                             fracice,     dqdt,      nsrflx,     qsrflx,  &
                             xx_mfup_max, xx_wcldbase, xx_kcldbase,       &
                             dcondt_resusp3d, conu2,  dcondt2,           &
                             state%ncol, pver,       ncnstaer,   nbins,   &
                             pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
                             convproc_do_evaprain_atonce,                 &
                             convproc_pom_spechygro,                      &
                             errmsg,     errflg )
      if (errflg /= 0) call endrun(trim(errmsg))

      call outfld( 'DP_MFUP_MAX', xx_mfup_max, pcols, lchnk )
      call outfld( 'DP_WCLDBASE', xx_wcldbase, pcols, lchnk )
      call outfld( 'DP_KCLDBASE', xx_kcldbase, pcols, lchnk )

      ! WETC = wet-deposition tendency, CONU = updraft mixing ratio (interstitial
      !    and cloud-borne); computed in aero_convproc_run and returned as out-args
      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)

            call outfld( trim(cnst_name_extd(1,mm))//'WETC', dcondt2(:,:,1,mm), pcols, lchnk )
            call outfld( trim(cnst_name_extd(1,mm))//'CONU', conu2(:,:,1,mm), pcols, lchnk )
            call outfld( trim(cnst_name_extd(2,mm))//'WETC', dcondt2(:,:,2,mm), pcols, lchnk )
            call outfld( trim(cnst_name_extd(2,mm))//'CONU', conu2(:,:,2,mm), pcols, lchnk )

         end do
      end do

   end subroutine aero_convproc_dp_intr

!=========================================================================================

end module aero_convproc_cam
