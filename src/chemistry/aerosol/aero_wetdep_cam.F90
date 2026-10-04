module aero_wetdep_cam
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use camsrfexch,    only: cam_out_t
  use physics_buffer,only: physics_buffer_desc, pbuf_get_index, pbuf_set_field, pbuf_get_field
  use constituents,  only: pcnst, cnst_name, cnst_get_ind
  use phys_control,  only: phys_getopts
  use ppgrid,        only: pcols, pver
  use physconst,     only: gravit

  use cam_abortutils,only: endrun
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc
  use infnan,        only: nan, assignment(=)

  use cam_history,   only: addfld, add_default, horiz_only, outfld
  use wetdep_cam,    only: wetdep_init

  use aerosol_properties_mod, only: aero_name_len
  use aerosol_properties_mod, only: aerosol_properties

  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use aerosol_instances_mod, only: aerosol_instances_get_state, &
                                   aerosol_instances_get_props, &
                                   aerosol_instances_get_num_models

  use aero_convproc_cam, only: aero_convproc_readnl, aero_convproc_init, aero_convproc_intr
  use aero_convproc_cam, only: convproc_do_evaprain_atonce
  use aero_convproc_cam, only: deepconv_wetdep_history

  use infnan,         only: nan, assignment(=)
  use perf_mod,       only: t_startf, t_stopf

  implicit none
  private

  public :: aero_wetdep_readnl
  public :: aero_wetdep_init
  public :: aero_wetdep_tend

  real(r8), parameter :: NOTSET = -huge(1._r8)
  real(r8) :: sol_facti_cloud_borne   = NOTSET
  real(r8) :: sol_factb_interstitial  = NOTSET
  real(r8) :: sol_factic_interstitial = NOTSET

  integer :: fracis_idx = -1
  integer :: rprddp_idx          = -1
  integer :: rprdsh_idx          = -1
  integer :: nevapr_shcu_idx     = -1
  integer :: nevapr_dpcu_idx     = -1

  logical :: wetdep_active = .false.
  integer :: nwetdep = 0
  logical :: convproc_do_aer = .false.
  logical,allocatable :: aero_cnst_lq(:,:)
  integer,allocatable :: aero_cnst_id(:,:)
  logical, public, protected :: wetdep_lq(pcnst) ! set flags true for constituents with non-zero tendencies

  integer :: nspec_max=0
  integer :: nele_tot            ! total number of aerosol elements
  class(aerosol_properties), pointer :: aero_props=>null()

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_wetdep_readnl(nlfile)
    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, masterprocid, mpi_character, mpi_real8, mpi_integer, mpi_success
    use spmd_utils,      only: mpi_logical

    character(len=*), intent(in)  :: nlfile  ! filepath for file containing namelist input

    integer                       :: unitn, ierr
    character(len=*), parameter   :: subname = 'aero_wetdep_readnl'

    ! ===================
    ! Namelist definition
    ! ===================
    namelist /aero_wetdep_nl/ sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial

    ! =============
    ! Read namelist
    ! =============
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aero_wetdep_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aero_wetdep_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)

       ! ============================
       ! Log namelist options
       ! ============================
       write(iulog,*) subname,' namelist settings: '
       write(iulog,*) '   sol_facti_cloud_borne  : ',sol_facti_cloud_borne
       write(iulog,*) '   sol_factb_interstitial : ',sol_factb_interstitial
       write(iulog,*) '   sol_factic_interstitial: ',sol_factic_interstitial
    end if

    ! ============================
    ! Broadcast namelist variables
    ! ============================
    call mpi_bcast(sol_facti_cloud_borne, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: sol_facti_cloud_borne')
    end if
    call mpi_bcast(sol_factb_interstitial, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: sol_factb_interstitial')
    end if
    call mpi_bcast(sol_factic_interstitial, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: sol_factic_interstitial')
    end if

    call mpi_bcast(nwetdep, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: nwetdep')
    end if

    wetdep_active = .true. !nwetdep>0

    if (masterproc) then
       write(iulog,*) subname,' wetdep_active = ',wetdep_active,' nwetdep = ',nwetdep
    endif

    call aero_convproc_readnl(nlfile)

  end subroutine aero_wetdep_readnl

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_wetdep_init( )

    use wetdep,       only: init_bcscavcoef
    use mo_constants, only: pi, boltz_cgs, rgas_cgs

    character(len=*), parameter :: subrname = 'aero_wetdep_init'

    character(len=2)  :: unit_basename  ! Units 'kg' or '1'
    character(len=aero_name_len) :: tmpname
    character(len=aero_name_len) :: tmpname_cw

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_chemistry

    integer :: l, m, id, astat, iaermod
    character(len=2) :: binstr
    class(aerosol_properties), pointer :: props_tmp

    character(len=512) :: errmsg   ! error handling for the portable init_bcscavcoef
    integer            :: errflg

    fracis_idx = pbuf_get_index('FRACIS')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    rprdsh_idx      = pbuf_get_index('RPRDSH')
    nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')

    if (.not.wetdep_active) return

    call phys_getopts(history_aerosol_out = history_aerosol, &
                      history_chemistry_out=history_chemistry, &
                      convproc_do_aer_out = convproc_do_aer)

    ! get the persistent properties for the aerosol model wetdep is working on
    ! (wetdep does not work with BAM)
    nullify(aero_props)
    do iaermod = 1, aerosol_instances_get_num_models()
       props_tmp => aerosol_instances_get_props(iaermod, 0)
       if (associated(props_tmp)) then
          if (.not. props_tmp%model_is('BAM')) then
             aero_props => props_tmp
             exit
          end if
       end if
    end do
    if (.not.associated(aero_props)) then
       call endrun(subrname//' : no non-BAM aerosol properties instance available')
    end if

    nele_tot = aero_props%ncnst_tot()

    allocate(aero_cnst_lq(aero_props%nbins(),0:maxval(aero_props%nspecies())), stat=astat)
    if (astat/=0) then
       call endrun(subrname//' : not able to allocate aero_cnst_lq array')
    end if
    aero_cnst_lq(:,:) = .false.

    allocate(aero_cnst_id(aero_props%nbins(),0:maxval(aero_props%nspecies())), stat=astat)
    if (astat/=0) then
       call endrun(subrname//' : not able to allocate aero_cnst_id array')
    end if
    aero_cnst_id(:,:) = -1

    wetdep_lq = .false.

    do m = 1, aero_props%nbins()
       write(binstr,'(i2.2)') m
       call addfld('SOLFACTB'//binstr,  (/ 'lev' /), 'A', '1', 'below cld sol fact')

       do l = 0, aero_props%nspecies(m)

          if (l == 0) then   ! number
             call aero_props%num_names( m, tmpname, tmpname_cw)
          else
             call aero_props%mmr_names( m,l, tmpname, tmpname_cw)
          end if

          call cnst_get_ind(tmpname, id, abort=.false.)
          aero_cnst_id(m,l) = id
          aero_cnst_lq(m,l) = id > 0
          if (id > 0) then
             wetdep_lq(id) = .true.
          end if

          ! units --
          if (l==0) then
             unit_basename = ' 1' ! for num
          else
             unit_basename = 'kg'
          endif

          call add_hist_fields(tmpname, unit_basename)
          call add_hist_fields(tmpname_cw, unit_basename)

          call addfld( trim(tmpname_cw)//'RSPTD', (/ 'lev' /), 'A', unit_basename//'/kg/s',   &
                       trim(tmpname_cw)//' resuspension tendency')

       end do
    end do

    call wetdep_init()

    nspec_max = maxval(aero_props%nspecies()) + 2

    ! build the below-cloud impaction/interception scavenging lookup table
    ! (allocation + fill are owned by the portable init_bcscavcoef)
    call init_bcscavcoef( aero_props, pi, boltz_cgs, rgas_cgs, &
                          errmsg, errflg )
    if (errflg /= 0) call endrun(trim(errmsg))

    if (convproc_do_aer) then
       call aero_convproc_init(aero_props)
    end if

  contains

    subroutine add_hist_fields(name,baseunits)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: baseunits

      call addfld (trim(name)//'SFWET', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux at surface')
      call addfld (trim(name)//'SFWETC', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (convective) at surface')
      call addfld (trim(name)//'SFSIC', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (incloud, convective) at surface')
      call addfld (trim(name)//'SFSIS', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (incloud, stratiform) at surface')
      call addfld (trim(name)//'SFSBC', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (belowcloud, convective) at surface')
      call addfld (trim(name)//'SFSBS', &
           horiz_only,  'A',baseunits//'/m2/s ','Wet deposition flux (belowcloud, stratiform) at surface')

      if (convproc_do_aer) then
         call addfld (trim(name)//'SFSEC', &
              horiz_only,  'A',unit_basename//'/m2/s','Wet deposition flux (precip evap, convective) at surface')
         call addfld (trim(name)//'SFSES', &
              horiz_only,  'A',unit_basename//'/m2/s','Wet deposition flux (precip evap, stratiform) at surface')
         call addfld (trim(name)//'SFSBD', &
              horiz_only,  'A',unit_basename//'/m2/s','Wet deposition flux (belowcloud, deep convective) at surface')
         call addfld (trim(name)//'WETC',  &
              (/ 'lev' /), 'A',unit_basename//'/kg/s ','wet deposition tendency')
         call addfld (trim(name)//'CONU',  &
              (/ 'lev' /), 'A',unit_basename//'/kg ','updraft mixing ratio')
      end if

      call addfld (trim(name)//'WET',(/ 'lev' /), 'A',baseunits//'/kg/s ','wet deposition tendency')
      call addfld (trim(name)//'INS',(/ 'lev' /), 'A',baseunits//'/kg/s ','insol frac')

      call addfld (trim(name)//'SIC',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' ic wet deposition')
      call addfld (trim(name)//'SIS',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' is wet deposition')
      call addfld (trim(name)//'SBC',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' bc wet deposition')
      call addfld (trim(name)//'SBS',(/ 'lev' /), 'A',baseunits//'/kg/s ', &
           trim(name)//' bs wet deposition')

      if ( history_aerosol .or. history_chemistry ) then
         call add_default (trim(name)//'SFWET', 1, ' ')
      endif
      if ( history_aerosol ) then
         call add_default (trim(name)//'SFSIC', 1, ' ')
         call add_default (trim(name)//'SFSIS', 1, ' ')
         call add_default (trim(name)//'SFSBC', 1, ' ')
         call add_default (trim(name)//'SFSBS', 1, ' ')
         if (convproc_do_aer) then
            call add_default (trim(name)//'SFSEC', 1, ' ')
            call add_default (trim(name)//'SFSES', 1, ' ')
            call add_default (trim(name)//'SFSBD', 1, ' ')
         end if
      endif

    end subroutine add_hist_fields

  end subroutine aero_wetdep_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_wetdep_tend( state, dt, dlf, cam_out, ptend, pbuf)
    use wetdep,     only: wetdepa_v2, get_bcscavcoefs
    use wetdep_cam, only: wetdep_inputs_set, wetdep_inputs_t
    use aerodep_flx, only: aerodep_flx_prescribed
    use aero_deposition_cam, only: aero_deposition_cam_setwet

    type(physics_state), target, intent(in) :: state  ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    character(len=*), parameter :: subrname = 'aero_wetdep_tend'
    type(wetdep_inputs_t) :: dep_inputs
    character(len=512) :: errmsg   ! error handling for the portable wetdepa_v2
    integer            :: errflg
    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble (pcols, pver, pcnst)
    real(r8), target :: fracis_nadv(pcols,pver)  ! fraction of not-transported aerosols

    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! interstitial num (1), interstitial vol (2)
    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: strt_loop, end_loop, stride_loop !loop indices for the lphase loop

    real(r8) :: sol_factb(pcols, pver)
    real(r8) :: sol_facti(pcols, pver)
    real(r8) :: sol_factic(pcols,pver)

    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species
    real(r8) :: rcscavt(pcols, pver)
    real(r8) :: rsscavt(pcols, pver)
    real(r8) :: iscavt(pcols, pver)
    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)

    real(r8) :: diam_wet(state%ncol, pver)
    logical  :: isprx(pcols,pver) ! true if precipation
    real(r8) :: prec(pcols) ! precipitation rate

    real(r8) :: rtscavt(pcols, pver, 0:nspec_max)

    integer :: ncol, lchnk, m, ndx,mm, l
    integer :: i,k

    real(r8), pointer :: insolfr_ptr(:,:)
    real(r8) :: q_tmp(pcols,pver) ! temporary array to hold "most current" mixing ratio for 1 species
    logical :: cldbrn

    type(ptr2d_t) :: raer(nele_tot)
    type(ptr2d_t) :: qqcw(nele_tot)

    real(r8) :: sflx(pcols)
    character(len=aero_name_len) :: aname, cname, name

    real(r8) :: qqcw_in(pcols,pver), qqcw_sav(pcols,pver,0:nspec_max)
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01

    character(len=2) :: binstr
    real(r8) :: aerdepwetcw(pcols,pcnst) ! aerosol wet deposition (cloud water)
    real(r8) :: aerdepwetis(pcols,pcnst)  ! aerosol wet deposition (interstitial)
    real(r8) :: dcondt_resusp3d(nele_tot,pcols,pver)

    integer, parameter :: nsrflx_mzaer2cnvpr = 2
    real(r8) :: qsrflx_mzaer2cnvpr(pcols,nele_tot,nsrflx_mzaer2cnvpr)

    real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
    real(r8), pointer :: rprdsh(:,:)     ! rain production, shallow convection
    real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.
    real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.

    real(r8) :: rprddpsum(pcols)
    real(r8) :: rprdshsum(pcols)
    real(r8) :: evapcdpsum(pcols)
    real(r8) :: evapcshsum(pcols)

    real(r8) :: tmp_resudp, tmp_resush
    real(r8) :: tmpa, tmpb
    real(r8) :: sflxec(pcols), sflxecdp(pcols)  ! deposition flux
    real(r8) :: sflxic(pcols), sflxicdp(pcols)  ! deposition flux
    real(r8) :: sflxbc(pcols), sflxbcdp(pcols)  ! deposition flux

    class(aerosol_state), pointer :: aero_state
    class(aerosol_properties), pointer :: props_tmp
    integer :: iaermod

    nullify(aero_state)

    if (.not.wetdep_active) return

    dcondt_resusp3d(:,:,:) = 0._r8

    !REMOVECAM - get persistent state from factory; under CAM-SIMA states will be passed as scheme inputs
    nullify(aero_state)
    do iaermod = 1, aerosol_instances_get_num_models()
       props_tmp => aerosol_instances_get_props(iaermod, 0)
       if (associated(props_tmp)) then
          if (.not. props_tmp%model_is('BAM')) then
             aero_state => aerosol_instances_get_state(iaermod, 0, state%lchnk)
             exit
          end if
       end if
    end do
    if (.not.associated(aero_state)) then
       call endrun(subrname//' : no non-BAM aerosol state available for wetdep')
    end if
    !REMOVECAM_END

    lchnk = state%lchnk
    ncol = state%ncol

    call physics_ptend_init(ptend, state%psetcols, subrname, lq=wetdep_lq)

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    call pbuf_get_field(pbuf, fracis_idx, fracis)

    call aero_state%get_states( aero_props, raer, qqcw )

    qsrflx_mzaer2cnvpr(:,:,:) = 0.0_r8
    aerdepwetis(:,:) = 0.0_r8
    aerdepwetcw(:,:) = 0.0_r8

    if (convproc_do_aer) then
       !Do cloudborne first for unified convection scheme so that the resuspension of cloudborne
       !can be saved then applied to interstitial
       strt_loop   =  2
       end_loop    =  1
       stride_loop = -1
    else
       ! Counters for "without" unified convective treatment (i.e. default case)
       strt_loop   = 1
       end_loop    = 2
       stride_loop = 1
    endif

    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
                    *state%pdel(:ncol,k)/gravit
    end do

    f_act_conv = 0._r8
    scavcoefnv = nan
    qqcw_sav = nan

    if (convproc_do_aer) then

       call t_startf('aero_convproc')
       call aero_convproc_intr( aero_props, aero_state, state, ptend, pbuf, dt, &
            nsrflx_mzaer2cnvpr, qsrflx_mzaer2cnvpr, aerdepwetis, dcondt_resusp3d )

       if (convproc_do_evaprain_atonce) then

          do m = 1,aero_props%nbins()
             do l = 0,aero_props%nspecies(m)
                mm = aero_props%indexer(m,l)

                if (l == 0) then   ! number
                   call aero_props%num_names(m, aname, cname)
                else
                   call aero_props%mmr_names(m,l, aname, cname)
                end if

                call outfld( trim(cname)//'RSPTD', dcondt_resusp3d(mm,:ncol,:), ncol, lchnk )

                do k = 1,pver
                   do i = 1,ncol
                      qqcw(mm)%fld(i,k) = max(0._r8, qqcw(mm)%fld(i,k) + dcondt_resusp3d(mm,i,k)*dt)
                   end do
                end do

             end do
          end do
       end if
       call t_stopf('aero_convproc')

    end if

    bins_loop: do m = 1,aero_props%nbins()

       phase_loop: do lphase = strt_loop,end_loop, stride_loop ! loop over interstitial (1) and cloud-borne (2) forms

          cldbrn = lphase==2

          sol_factb = nan
          sol_facti = nan
          sol_factic = nan

          if (lphase == 1) then ! interstial aerosol

             sol_facti = 0.0_r8 ! strat in-cloud scav totally OFF for institial

             sol_factic = sol_factic_interstitial

          else ! cloud-borne aerosol (borne by stratiform cloud drops)

             sol_factb  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti  = sol_facti_cloud_borne   ! strat  in-cloud scav cloud-borne tuning factor
             sol_factic = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
             !        that conv precip collects strat droplets)
             f_act_conv = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean

          end if
          if (convproc_do_aer .and. lphase == 1) then
             ! if modal aero convproc is turned on for aerosols, then
             !    turn off the convective in-cloud removal for interstitial aerosols
             !    (but leave the below-cloud on, as convproc only does in-cloud)
             !    and turn off the outfld SFWET, SFSIC, SFSID, SFSEC, and SFSED calls
             ! for (stratiform)-cloudborne aerosols, convective wet removal
             !    (all forms) is zero, so no action is needed
             sol_factic = 0.0_r8
          endif

          diam_wet = aero_state%wet_diameter(m,ncol,pver)

          scavcoefnv = 0.0_r8

          if (lphase == 1) then ! interstial aerosol
             call get_bcscavcoefs( m, ncol, pver, isprx, diam_wet, scavcoefnv(:,:,1), scavcoefnv(:,:,2), aero_props )

             if ( sol_factb_interstitial /= NOTSET ) then
                sol_factb(:ncol,:) = sol_factb_interstitial ! all below-cloud scav
             else
                sol_factb(:ncol,:) = aero_state%sol_factb_interstitial( m, ncol, pver, aero_props )
             end if

             write(binstr,'(i2.2)') m
             call outfld('SOLFACTB'//binstr,sol_factb, pcols, lchnk)

          end if

          elem_loop: do l = 0,aero_props%nspecies(m)

             ndx = aero_cnst_id(m,l)
             if (ndx<1) cycle elem_loop

             if (.not. cldbrn .and. ndx>0) then
                insolfr_ptr => fracis(:,:,ndx)
             else
                insolfr_ptr => fracis_nadv
             endif

             mm = aero_props%indexer(m,l)

             if (l == 0) then   ! number
                call aero_props%num_names( m, aname, cname)
             else
                call aero_props%mmr_names( m,l, aname, cname)
             end if

             if (cldbrn) then
                q_tmp(1:ncol,:) = qqcw(mm)%fld(1:ncol,:)
                jnv = 0
                if (convproc_do_aer) then
                   qqcw_sav(:ncol,:,l) = q_tmp(1:ncol,:)
                endif
                name = cname
                qqcw_in = nan
                f_act_conv = nan
             else ! interstial aerosol
                q_tmp(1:ncol,:) = raer(mm)%fld(1:ncol,:) + ptend%q(1:ncol,:,ndx)*dt
                if (l==0) then
                   jnv = 1
                else
                   jnv = 2
                end if
                if(convproc_do_aer) then
                   !Feed in the saved cloudborne mixing ratios from phase 2
                   qqcw_in(:ncol,:) = qqcw_sav(:ncol,:,l)
                else
                   qqcw_in(:ncol,:) = qqcw(mm)%fld(:ncol,:)
                end if

                f_act_conv(:ncol,:) = aero_state%convcld_actfrac( aero_props, m, l, ncol, pver)
                name = aname
             end if

             dqdt_tmp(1:ncol,:) = 0.0_r8

             call wetdepa_v2(state%pdel, &
                  dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                  dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, &
                  dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                  dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                  dlf, insolfr_ptr, sol_factb(:ncol,:), ncol, &
                  scavcoefnv(:,:,jnv), gravit, pver, errmsg, errflg, &
                  is_strat_cloudborne=cldbrn, &
                  qqcw=qqcw_in(:,:), f_act_conv=f_act_conv, &
                  icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                  convproc_do_aer=convproc_do_aer, rcscavt=rcscavt, rsscavt=rsscavt,  &
                  sol_facti_in=sol_facti(:ncol,:), sol_factic_in=sol_factic(:ncol,:), &
                  convproc_do_evaprain_atonce_in=convproc_do_evaprain_atonce, &
                  bergso_in=dep_inputs%bergso )
             if (errflg /= 0) call endrun(trim(errmsg))

             if(convproc_do_aer) then
                if(cldbrn) then
                   ! save resuspension of cloudborne species
                   rtscavt(1:ncol,:,l) = rcscavt(1:ncol,:) + rsscavt(1:ncol,:)
                   ! wetdepa_v2 adds the resuspension of cloudborne to the dqdt of cloudborne (as a source)
                   ! undo this, so the resuspension of cloudborne can be added to the dqdt of interstitial (above)
                   dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) - rtscavt(1:ncol,:,l)
                else
                   ! add resuspension of cloudborne species to dqdt of interstitial species
                   dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) + rtscavt(1:ncol,:,l)
                end if
             endif

             if (cldbrn) then
                do k = 1,pver
                   do i = 1,ncol
                      if ( (qqcw(mm)%fld(i,k) + dqdt_tmp(i,k) * dt) .lt. 0.0_r8 )   then
                         dqdt_tmp(i,k) = - qqcw(mm)%fld(i,k) / dt
                      end if
                   end do
                end do

                qqcw(mm)%fld(1:ncol,:) = qqcw(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

             else
                ptend%q(1:ncol,:,ndx) = ptend%q(1:ncol,:,ndx) + dqdt_tmp(1:ncol,:)
             end if

             call outfld( trim(name)//'WET', dqdt_tmp(:,:), pcols, lchnk)
             call outfld( trim(name)//'SIC', icscavt, pcols, lchnk)
             call outfld( trim(name)//'SIS', isscavt, pcols, lchnk)
             call outfld( trim(name)//'SBC', bcscavt, pcols, lchnk)
             call outfld( trim(name)//'SBS', bsscavt, pcols, lchnk)

             call outfld( trim(name)//'INS', insolfr_ptr, pcols, lchnk)

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                enddo
             enddo
             if (cldbrn) then
                call outfld( trim(name)//'SFWET', sflx, pcols, lchnk)
                if (ndx>0) aerdepwetcw(:ncol,ndx) = sflx(:ncol)
             else
                if (ndx>0) aerdepwetis(:ncol,ndx) = aerdepwetis(:ncol,ndx) + sflx(:ncol)
                call outfld( trim(name)//'SFWET', aerdepwetis(:ncol,ndx), ncol, lchnk)
             end if

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                enddo
             enddo
             if (cldbrn) then
                call outfld( trim(name)//'SFSIC', sflx, pcols, lchnk)
             else
                if (.not.convproc_do_aer) call outfld( trim(name)//'SFSIC', sflx, pcols, lchnk)
                if (convproc_do_aer) sflxic = sflx
             end if

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(name)//'SFSIS', sflx, pcols, lchnk)

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(name)//'SFSBC', sflx, pcols, lchnk)
             if (convproc_do_aer) sflxbc = sflx

             sflx(:)=0._r8
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                enddo
             enddo
             call outfld( trim(name)//'SFSBS', sflx, pcols, lchnk)

             if(convproc_do_aer) then

                sflx(:)=0.0_r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+rsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(name)//'SFSES', sflx, pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+rcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                if (.not.convproc_do_aer) call outfld( trim(name)//'SFSEC', sflx, pcols, lchnk)
                sflxec = sflx

                if(.not.cldbrn) then

                   ! apportion convective surface fluxes to deep and shallow conv
                   ! this could be done more accurately in subr wetdepa
                   ! since deep and shallow rarely occur simultaneously, and these
                   !    fields are just diagnostics, this approximate method is adequate
                   ! only do this for interstitial aerosol, because conv clouds to not
                   !    affect the stratiform-cloudborne aerosol
                   if ( deepconv_wetdep_history) then

                      call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
                      call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
                      call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
                      call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )

                      rprddpsum(:)  = 0.0_r8
                      rprdshsum(:)  = 0.0_r8
                      evapcdpsum(:) = 0.0_r8
                      evapcshsum(:) = 0.0_r8

                      do k = 1, pver
                         rprddpsum(:ncol)  = rprddpsum(:ncol)  +  rprddp(:ncol,k)*state%pdel(:ncol,k)/gravit
                         rprdshsum(:ncol)  = rprdshsum(:ncol)  +  rprdsh(:ncol,k)*state%pdel(:ncol,k)/gravit
                         evapcdpsum(:ncol) = evapcdpsum(:ncol) + evapcdp(:ncol,k)*state%pdel(:ncol,k)/gravit
                         evapcshsum(:ncol) = evapcshsum(:ncol) + evapcsh(:ncol,k)*state%pdel(:ncol,k)/gravit
                      end do

                      do i = 1, ncol
                         rprddpsum(i)  = max( rprddpsum(i),  1.0e-35_r8 )
                         rprdshsum(i)  = max( rprdshsum(i),  1.0e-35_r8 )
                         evapcdpsum(i) = max( evapcdpsum(i), 0.1e-35_r8 )
                         evapcshsum(i) = max( evapcshsum(i), 0.1e-35_r8 )

                         ! assume that in- and below-cloud removal are proportional to column precip production
                         tmpa = rprddpsum(i) / (rprddpsum(i) + rprdshsum(i))
                         tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )
                         sflxicdp(i) = sflxic(i)*tmpa
                         sflxbcdp(i) = sflxbc(i)*tmpa

                         ! assume that resuspension is proportional to (wet removal)*[(precip evap)/(precip production)]
                         tmp_resudp =           tmpa  * min( (evapcdpsum(i)/rprddpsum(i)), 1.0_r8 )
                         tmp_resush = (1.0_r8 - tmpa) * min( (evapcshsum(i)/rprdshsum(i)), 1.0_r8 )
                         tmpb = max( tmp_resudp, 1.0e-35_r8 ) / max( (tmp_resudp+tmp_resush), 1.0e-35_r8 )
                         tmpb = max( 0.0_r8, min( 1.0_r8, tmpb ) )
                         sflxecdp(i) = sflxec(i)*tmpb
                      end do
                      call outfld( trim(name)//'SFSBD', sflxbcdp, pcols, lchnk)
                   else
                      sflxec(1:ncol)   = 0.0_r8
                      sflxecdp(1:ncol) = 0.0_r8
                   end if

                   ! when ma_convproc_intr is used, convective in-cloud wet removal is done there
                   ! the convective (total and deep) precip-evap-resuspension includes in- and below-cloud
                   ! contributions
                   ! so pass the below-cloud contribution to ma_convproc_intr
                   qsrflx_mzaer2cnvpr(1:ncol,mm,1) = sflxec(  1:ncol)
                   qsrflx_mzaer2cnvpr(1:ncol,mm,2) = sflxecdp(1:ncol)

                end if
             end if

          end do elem_loop
       end do phase_loop

    end do bins_loop

    nullify(aero_state)

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not. aerodep_flx_prescribed()) then
       call aero_deposition_cam_setwet(aerdepwetis, aerdepwetcw, cam_out)
    endif

  end subroutine aero_wetdep_tend

end module aero_wetdep_cam
